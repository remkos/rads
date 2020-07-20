!-----------------------------------------------------------------------
! Copyright (c) 2011-2020  Remko Scharroo
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

!*radspassesindex -- Print information about passes
!+
program radspassesindex

! This is a little helper program that prints out the following
! information per pass:
! - Phase name
! - Cycle number
! - Pass number
! - Start time
! - End time
! - Equator crossing time
! - Equator crossing longitude
! - Number of points
!
! Use -h to get a header, -p to predict the values (without reading)
!-----------------------------------------------------------------------
use rads
use rads_time

integer(fourbyteint) :: i, cycle, pass
logical :: header = .false., predict = .false., ymd = .false.
character(len=26) :: date(3) = &
(/ 'start_time                ', 'end_time                  ', 'equator_time              ' /)
type(rads_sat) :: S
type(rads_pass) :: P

! Initialise RADS or issue help
call synopsis
call rads_set_options ('hpy')
call rads_init (S)
if (S%error /= rads_noerr) call rads_exit ('Fatal error')

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('h')
		header = .true.
	case ('p')
		predict = .true.
	case ('y')
		ymd = .true.
	end select
enddo

! Print header if requested
if (.not.header) then
	! Skip
else if (ymd) then
	write (*,600) date
else
	write (*,600) date(:)(1:17)
endif

! Estimate the equator crossings
do cycle = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		if (predict) then
			call rads_predict_equator (S, P, cycle, pass)
			if (S%error /= rads_noerr) cycle
			P%ndata = floor((P%end_time - P%start_time)/S%dt1hz + 1)
		else
			call rads_open_pass (S, P, cycle, pass)
		endif
		P%equator_lon = modulo(P%equator_lon, 360d0)
		if (P%ndata <= 0) then
			! Skip
		else if (ymd) then
			date(1) = strf1985f(P%start_time, 'T')
			date(2) = strf1985f(P%end_time, 'T')
			date(3) = strf1985f(P%equator_time, 'T')
			write (*,611) S%phase%name, cycle, pass, date, P%equator_lon, P%ndata
		else
			write (*,610) S%phase%name, cycle, pass, P%start_time, P%end_time, P%equator_time, P%equator_lon, P%ndata
		endif
		if (.not.predict) call rads_close_pass (S, P, .true.)
	enddo
enddo

call rads_end (S)

! Formats
600 format ('# cyc pass ',3(a,1x),'equator_lon  nr')
610 format (a1,1x,i3.3,1x,i4.4,3f18.6,f11.6,i5)
611 format (a1,1x,i3.3,1x,i4.4,3(1x,a26),f11.6,i5)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('Print information about passes')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -h                        Add header'/ &
'  -p                        Predict the pass information rather than opening files'/ &
'  -y                        Print out dates as YYYY-MM-DDTHH:MM:SS.SSSSSS')
stop
end subroutine synopsis

!***********************************************************************

end program radspassesindex
