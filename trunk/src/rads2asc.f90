!-----------------------------------------------------------------------
! $Id$
!
! Copyright (C) 2012  Remko Scharroo (Altimetrics LLC)
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
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
integer(fourbyteint), parameter :: msat = 5
type(rads_sat), target :: Sats(msat)
type(rads_sat), pointer :: S
type(rads_pass) :: P
type(rads_var), pointer :: var, temp(:)

! Local declarations, etc.
integer(fourbyteint) :: outunit, listunit = -1, logunit = 6
character(len=rads_naml) :: outname = ''
character(len=640) :: format_string
integer(fourbyteint) :: i, j, l, ios, cycle, pass, step = 1, nseltot = 0, nselmax = huge(0_fourbyteint)
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
call rads_set_options ('r::fs:o: step: list: output: maxrec:')
call rads_init (Sats)
if (any(Sats%error /= rads_noerr)) call rads_exit ('Fatal error')

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('o', 'output')
		outname = rads_opt(i)%arg
		if (outname == '') outname = '-'
	case ('r')
		if (rads_opt(i)%arg == 'n') then
			reject = -2
		else
			reject = 0
			read (rads_opt(i)%arg, *, iostat=ios) reject
		endif
	case ('f')
		freeform = .true.
	case ('s')
		if (rads_opt(i)%arg == 'p') then
			stat_mode = pass_stat
		else if (rads_opt(i)%arg == 'c') then
			stat_mode = cycle_stat
		else
			stat_mode = no_stat
		endif
	case ('maxrec')
		read (rads_opt(i)%arg, *, iostat=ios) nselmax
	case ('step')
		read (rads_opt(i)%arg, *, iostat=ios) step
	case ('list')
		listunit = getlun()
		open (listunit, file=rads_opt(i)%arg, status='replace')
	end select
enddo

! Now loop through all satellites

do j = 1,msat
	if (Sats(j)%sat == '') exit
	S => Sats(j)

	! If there is no -f option, add relative time, lat, lon to the front of the list
	if (.not.freeform) then
		allocate (temp(S%nsel+3))
		var => rads_varptr (S, 'time_rel_eq')
		temp(1) = var
		temp(2) = S%lat
		temp(3) = S%lon
		temp(4:S%nsel+3) = S%sel(1:S%nsel)
		S%nsel = S%nsel + 3
		deallocate (S%sel, stat=ios)
		S%sel => temp
	endif

	! Per-cycle statistics cannot be done for per-pass files
	if (outname == '' .and. stat_mode == cycle_stat) then
		call rads_error (S, rads_noerr, 'Cannot do per-cycle statistics for per-pass files. -sc ignored')
		stat_mode = no_stat
	endif

	! If flags or SLA are among the results, remember which they are
	do i = 1,S%nsel
		if (S%sel(i)%info%datatype == rads_type_flagmasks) then
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
760 format(/'Maximum number of output records reached (',i0,' >= ',i0,')')

	! Finish progress bar
	if (S%debug >= 1) write (logunit,*)

	! Close file before exit
	if (outunit /= 6) close (outunit)

	! Deallocate temporary arrays
	deallocate (stat, r, q)
	! No need to deallocate temp ... will be done by rads_end.

enddo	! Next satellite

! Print overall statistics and close RADS
call rads_stat (Sats, logunit)
call rads_end (Sats)

contains

!***********************************************************************

subroutine synopsis ()
if (rads_version ('$Revision$','Select RADS altimeter data and output to ASCII')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r#                       Reject lines if data item number # on sel= specifier is NaN'/ &
'                            (default: reject if SLA field is NaN)'/ &
'  -r0, -r                   Do not reject lines with NaN values'/ &
'  -rn                       Reject lines if any value is NaN'/ &
'  -f                        Do not start with t,lat,lon in output'/ &
'  -sp, -sc                  Include statistics per pass or per cycle'/ &
'  --step=N                  Step through records with stride N (default = 1)'/ &
'  --list=LISTNAME           Specify creation of list of output files'/ &
'  --maxrec=N                Specify maximum number of output records (default = unlimited)'/ &
'  -o, --out=OUTNAME         Specify name of a single output file (default is pass files, - is stdout)')
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
'# Created: ',a,' UTC: ',a/&
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
write (outunit,600) timestamp(), trim(S%command), trim(S%satellite), trim(S%phase%name), cycle, pass, &
	P%equator_time, strf1985f (P%equator_time), P%equator_lon, trim(P%original)

! Write column info
if (outname /= '') continued = .true.
m = 0
do j = 1,S%nsel
	if (S%sel(j)%info%datatype == rads_type_flagmasks) then
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
