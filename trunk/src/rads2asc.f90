!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2011  Remko Scharroo (Altimetrics LLC)
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!-----------------------------------------------------------------------

!*rads2asc -- Select RADS data and output to ASCII
!+
program rads2asc
!
! This program converts the RADS netCDF altimeter data to ASCII.
! At the same time it applies the standard selection criteria
! and allows some further modifications and selections.
!
! usage: rads2asc sat=<sat> [RADS_options] [options]
!-----------------------------------------------------------------------
use rads
use rads_netcdf
use rads_misc
use rads_time

! RADS structures
type(rads_sat) :: S
type(rads_pass) :: P
type(rads_var), pointer :: temp(:)

! Local declarations, etc.
integer(fourbyteint) :: outunit, listunit = -1, logunit = 6
character(len=256) :: arg, outname = ''
character(len=512) :: format_string
integer(fourbyteint) :: i, l, ios, cycle, pass, step = 1, nseltot = 0, nselmax = huge(0_fourbyteint)
logical :: freeform = .false.
integer(fourbyteint) :: isflags = 0, reject = -1

! For statistics
integer, parameter :: no_stat = 0, pass_stat = 1, cycle_stat = 2
integer :: stat_mode = no_stat, nselpass = 0
type :: stat_type
	integer(fourbyteint) :: nr
	real(eightbytereal) :: mean, sum2, xmin, xmax
endtype
type(stat_type), allocatable :: stat(:)
real(eightbytereal), allocatable :: r(:), q(:)

! Initialize RADS or issue help
call synopsis
call rads_init (S)
if (S%error /= rads_noerr) call rads_exit ('Fatal error')

! Scan command line arguments
do i = 1,iargc()
	call getarg(i,arg)
	if (arg(:4) == 'out=') then
		outname = arg(5:)
	else if (arg(:2) == '-rn') then
		reject = -2
	else if (arg(:2) == '-r') then
		reject = 0
		read (arg(3:),*,iostat=ios) reject
	else if (arg(:2) == '-f') then
		freeform = .true.
	else if (arg(:3) == '-sp') then
		stat_mode = pass_stat
	else if (arg(:3) == '-sc') then
		stat_mode = cycle_stat
	else if (arg(:7) == 'maxrec=') then
		read (arg(8:),*) nselmax
	else if (arg(:5) == 'step=') then
		read (arg(6:),*) step
	else if (arg(:5) == 'list=') then
		listunit = getlun()
		open (listunit, file=arg(6:), status='replace')
	endif
enddo

! If there is no -f option, add relative time, lat, lon to the front of the list
if (.not.freeform) then
	allocate (temp(S%nsel+3))
	temp(4:S%nsel+3) = S%sel(1:S%nsel)
	deallocate (S%sel, stat=ios)
	S%sel => temp
	temp(1) = rads_varptr (S, 'time_rel_eq')
	temp(2) = rads_varptr (S, 'lat')
	temp(3) = rads_varptr (S, 'lon')
	S%nsel = S%nsel + 3
endif

! Per-cycle statistics cannot be done for per-pass files
if (outname == '' .and. stat_mode == cycle_stat) then
	call rads_error (S, rads_noerr, 'Cannot do per-cycle statistics for per-pass files. -sc ignored')
	stat_mode = no_stat
endif

! If flags or SLA are among the results, remember which they are
do i = 1,S%nsel
	if (S%sel(i)%info%datatype == rads_type_flagword) then
		isflags = i
	else if (S%sel(i)%info%datatype == rads_type_sla) then
		if (reject == -1) reject = i
	endif
enddo

! Open single output file, if requested
if (outname == '-') then
	outunit = 6
	logunit = 0
else if (outname /= '') then
	outunit = getlun()
	open (outunit, file=outname, status='replace')
	if (listunit >= 0) write (listunit,'(a)') trim(outname)
endif

! Determine format string for ASCII output
format_string = '('
l = 1
do i = 1,S%nsel
	format_string(l+1:) = trim(S%sel(i)%info%format) // ',1x,'
	l = len_trim (format_string)
enddo
format_string(l-3:) = ')'

! Initialise statistics
allocate (stat(S%nsel), r(S%nsel), q(S%nsel))
stat = stat_type (0, 0d0, 0d0, S%nan, S%nan)

! Now loop through all cycles and passes
do cycle = S%cycles(1), S%cycles(2), S%cycles(3)
	! Stop processing after too many output lines
	if (nseltot >= nselmax) then
		write (logunit,760) nseltot,nselmax
		exit
	endif

	! Process passes one-by-one
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cycle, pass)
		if (P%ndata > 0) call process_pass
		if (S%debug >= 1 .and. outunit /= logunit) call rads_progress_bar (S, P, nselpass, logunit)
		call rads_close_pass (S, P)
	enddo

	! Print out per-cycle statistics, it requested
	if (stat_mode == cycle_stat) call print_stat
enddo
760 format(/'Maximum number of output records reached (',i9,' >=',i9,')')

! Finish progress bar
if (S%debug >= 1) write (logunit,*)

! Close file before exit
if (outunit /= 6) close (outunit)

! Print overall statistics and close RADS
call rads_stat (S, logunit)
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('$Revision$','Select RADS altimeter data and output to ASCII')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r#               : reject lines if data item number # on sel= specifier is NaN'/ &
'                      (default: reject if SLA field is NaN)'/ &
'  -r0, -r           : do not reject lines with NaN values'/ &
'  -rn               : reject lines if any value is NaN'/ &
'  -f                : do not start with t,lat,lon in output'/ &
'  -sp, -sc          : include statistics per pass or per cycle'/ &
'  step=n            : step through records with stride n (default = 1)'/ &
'  list=listname     : specify creation of list of output files'/ &
'  out=outname       : specify name of a single output file (default is pass files, - is stdout)')
stop
end subroutine synopsis

!***********************************************************************
! Process data for a single pass

subroutine process_pass
character(len=80) :: passname
real(eightbytereal), allocatable :: data(:,:)
integer(fourbyteint) :: i

! Open a new output file when outname is not specified
if (outname == '') then
	write (passname, '(a2,"p",i4.4,"c",i3.3,".asc")') S%sat, pass, cycle
	outunit = getlun()
	open (outunit, file=passname, status='replace')
endif

! Read the data
allocate (data(P%ndata,S%nsel))

do i = 1,S%nsel
	call rads_get_var (S, P, S%sel(i), data(:,i))
enddo
nselpass = 0

! Loop through the data
do i = 1,P%ndata,step
	! See if we have to reject this record
	if (reject > 0) then
		if (isnan(data(i,reject))) cycle
	else if (reject == -2) then
		if (any(isnan(data(i,:)))) cycle
	endif
	nselpass = nselpass + 1

	! Update the pass-by-pass or cycle-by-cycle statistics
	if (stat_mode /= no_stat) then
		stat%nr = stat%nr + 1
		stat%xmin = min(stat%xmin, data(i,:))
		stat%xmax = max(stat%xmax, data(i,:))
		q = data(i,:) - stat%mean
		r = q / stat%nr
		stat%mean = stat%mean + r
		stat%sum2 = stat%sum2 + r * q * (stat%nr - 1)
	endif

	! Write header before first record of pass.
	! This placement of write_header makes sure the header is only written when there
	! are actual measurements passing the editing criteria.
	if (nselpass == 1) call write_header

	! Write out the data point, and take care of the flag words
	if (isflags > 0) then
		write (outunit,format_string) data(i,:isflags-1),to16bits(data(i,isflags)),data(i,isflags+1:)
	else if (S%nsel > 0) then
		write (outunit,format_string) data(i,:)
	endif
enddo

nseltot = nseltot + nselpass
deallocate (data)

! Print out pass statistics and reset them, if requested
if (stat_mode == pass_stat) call print_stat

! Close per-pass output file
if (outname /= '') then
	! Keep file open
else if (nselpass == 0) then
	close (outunit,status='delete')
else
	close (outunit,status='keep')
	if (listunit >= 0) write (listunit,'(a)') trim(passname)
endif

end subroutine process_pass

!***********************************************************************
! Print statistics for one batch of data

subroutine print_stat
if (stat(1)%nr == 0) return
call stat_line ('min ', stat%xmin)
call stat_line ('max ', stat%xmax)
call stat_line ('mean', stat%mean)
call stat_line ('std ', sqrt(stat%sum2/stat(1)%nr))
stat = stat_type (0, 0d0, 0d0, S%nan, S%nan)
end subroutine print_stat

!***********************************************************************

subroutine stat_line (string,x)
character(len=*) :: string
real(eightbytereal) :: x(:)
write (outunit, 690, advance='no') string, cycle, pass, stat(1)%nr
if (isflags > 0) then
	write (outunit, format_string) x(:isflags-1),to16bits(0d0),x(isflags+1:)
else
	write (outunit, format_string) x
endif
690 format('# ',a,t7,': ',i3.3,1x,i4.4,i9,' ')
end subroutine stat_line

!***********************************************************************

subroutine write_header
logical :: continued = .false.
integer :: j, k, m
character(len=4), parameter :: flagname(0:15) = (/'2516','2501','2502','2503','2504','2505','2506','2507', &
	'2508','2509','2510','2511','2512','2513','2514','2515'/)
type(rads_var), pointer :: var

! Format of ASCII header
600 format( &
'# RADS_ASC'/ &
'# Satellite = ',a/ &
'# Phase     = ',a/ &
'# Cycle     = ',i3.3/ &
'# Pass      = ',i4.4/ &
'# Equ_time  = ',f17.6,1x,'(',a,')'/ &
'# Equ_lon   = ',f11.6/ &
'# Original  = ',a)
620 format('# Col ',i2,'    = ',a,' [',a,']')
621 format('# Col ',i2,'    = flag ',i2,': ',a)

! Print the top of the header (skip a line first if this is a continuation)
if (continued) write (outunit,*)
write (outunit,600) trim(S%satellite), trim(S%phase%name), cycle, pass, &
	P%equator_time, strf1985f (P%equator_time), P%equator_lon, trim(P%original)

! Write column info
if (outname /= '') continued = .true.
m = 0
do j = 1,S%nsel
	if (S%sel(j)%info%datatype == rads_type_flagword) then
		do k = 0,15
			m = m + 1
			var => rads_varptr (S, flagname(k))
			write (outunit,621) m, k, trim(var%info%long_name)
		enddo
	else
		m = m + 1
		write (outunit,620) m, trim(S%sel(j)%info%long_name), trim(S%sel(j)%info%units)
	endif
enddo
end subroutine write_header

!***********************************************************************

end program rads2asc
