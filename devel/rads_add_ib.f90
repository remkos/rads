!-----------------------------------------------------------------------
! Copyright (c) 2011-2018  Remko Scharroo
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

!*rads_add_ib -- Add global inverse barometer correction to RADS data
!+
! This program adds the global inverse barometer correction field
! to the contents of RADS altimeter data files
!
! Note: The values for the global inverse barometer correction are
! based on tables of the mean global pressure over oceans based on
! the ECMWF or NCEP sea surface pressure grids. The source of the
! global mean values should match the models used for the inverse
! barometer correction. This generally is NCEP for GFO and GEOSAT
! and ECMWF for the other satellites.
!
! usage: rads_add_ib [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_ib

use rads
use rads_grid
use rads_misc
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Command line arguments

integer(fourbyteint) :: cyc, pass, type

! Initialise

call synopsis ('--head')
call rads_init (S)

! For GFO and GEOSAT use smoothed NCEP values. Otherwise use smoothed
! ECMWF values.

if (S%sat == 'g1' .or. S%sat == 'gs') then
	type = 4
else
	type = 2
endif

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

! Cleanup

call rads_end (S)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add global average inverse barometer correction to RADS data', flag=flag)) return
call synopsis_devel ('')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n), z(n)
integer(fourbyteint) :: i

call log_pass (P)

! Get time

call rads_get_var (S, P, 'time', time, .true.)

! Process all data points

do i = 1,n

! Get mean pressure at present time

	call globpres (type, time(i), z(i))
	if (z(i) == 1013.3d0) call globpres (6-type, time(i), z(i))

! Compute inverse barometer correction from global mean pressure.

	z(i) = -9.948d-3 * (1013.3d0-z(i))
enddo

! Store all data fields.

call rads_put_history (S, P)
call rads_def_var (S, P, 'inv_bar_global')
call rads_put_var (S, P, 'inv_bar_global', z)

call log_records (n)
end subroutine process_pass

end program rads_add_ib
