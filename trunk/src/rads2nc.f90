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

!*rads2nc -- Select RADS data and output to netCDF
!+
program rads2nc
!
! This program converts the RADS netCDF altimeter data and output them
! to new netCDF data files.
! At the same time it applies the standard selection criteria
! and allows some further modifications and selections.
!
! usage: rads2nc sat=<sat> [RADS_options] [options]
!-----------------------------------------------------------------------
use rads
use rads_netcdf
use rads_misc
use netcdf

! RADS structures
type(rads_sat) :: S
type(rads_pass) :: P, Pout

! Local declarations, etc.
integer(fourbyteint) :: logunit = 6
character(len=rads_naml) :: arg, opt, optarg, outname = ''
integer(fourbyteint) :: i, ios, cycle, pass, step = 1, nseltot = 0, nselmax = huge(0_fourbyteint)
integer(fourbyteint) :: reject = -1

! Output file definitions
integer(fourbyteint) :: nselpass = 0

! Initialize RADS or issue help
call synopsis
call rads_init (S)
if (S%error /= rads_noerr) call rads_exit ('Fatal error')

! Scan command line arguments
do i = 1,iargc()
	call getarg (i, arg)
	call splitarg (arg, opt, optarg)
	select case (opt)
	case ('-o', '--out')
		outname = optarg
	case ('-r')
		if (optarg == 'n') then
			reject = -2
		else
			reject = 0
			read (optarg, *, iostat=ios) reject
		endif
	case ('--maxrec')
		read (optarg, *, iostat=ios) nselmax
	case ('--step')
		read (optarg, *, iostat=ios) step
	end select
enddo

! If SLA is among the results, remember which index that is
do i = 1,S%nsel
	if (S%sel(i)%info%datatype == rads_type_sla) then
		if (reject == -1) reject = i
	endif
enddo

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
		if (S%debug >= 1) call rads_progress_bar (S, P, nselpass, logunit)
		call rads_close_pass (S, P)
	enddo
enddo
760 format(/'Maximum number of output records reached (',i0,' >= ',i0,')')

! Finish progress bar
if (S%debug >= 1) write (logunit,*)

! Close data file before exit
if (outname /= '') call rads_close_pass (S, Pout)

! Print overall statistics and close RADS
call rads_stat (S, logunit)
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('$Revsion 4.0 $','Select RADS altimeter data and output to netCDF')) return
call rads_synopsis ()
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r#                       Reject lines if data item number # on sel= specifier is NaN'/ &
'                            (default: reject if SLA field is NaN)'/ &
'  -r0, -r                   Do not reject lines with NaN values'/ &
'  -rn                       Reject lines if any value is NaN'/ &
'  --step=N                  Step through records with stride n (default = 1)'/ &
'  -o, --out=OUTNAME         Specify name of a single output file (default is pass files)')
stop
end subroutine synopsis

!***********************************************************************

subroutine process_pass
real(eightbytereal), allocatable :: data(:,:)
logical, allocatable :: accept(:)
integer(fourbyteint) :: i, start(1)

! Read the data
nselpass = 0
allocate (data(P%ndata,S%nsel), accept(P%ndata))
accept = .false.
do i = 1,S%nsel
	call rads_get_var (S, P, S%sel(i), data(:,i))
enddo

! Loop through the data
do i = 1,P%ndata,step
	! See if we have to reject this record
	if (reject > 0) then
		if (isnan(data(i,reject))) cycle
	else if (reject == -2) then
		if (any(isnan(data(i,:)))) cycle
	endif
	nselpass = nselpass + 1
	accept(i) = .true.
enddo

! If no data left, then return
if (nselpass == 0) return

! Open output pass file
start(1) = nseltot + 1
if (outname == '') then
	Pout = P
	call rads_create_pass (S, Pout, nselpass, '')
	call rads_def_var (S, Pout, S%sel)
	start(1) = 1
else if (nseltot == 0) then
	call rads_create_pass (S, Pout, 0, outname)
	call rads_def_var (S, Pout, S%sel)
endif

! Write the variables
do i = 1,S%nsel
	call rads_put_var (S, Pout, S%sel(i), pack(data(:,i),accept), start)
enddo
nseltot = nseltot + nselpass
deallocate (data, accept)

! Close per-pass output file
if (outname == '') call rads_close_pass (S, Pout)

end subroutine process_pass

!***********************************************************************

end program rads2nc
