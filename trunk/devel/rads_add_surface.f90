!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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

!*rads_add_surface -- Add surface type flags to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! based on a grid of surface type values from GMT/GSHHS.
! The program alters BOTH the flag word and the surface_type variable.
!
! usage: rads_add_surface [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_surface

use rads
use rads_grid
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
type(grid) :: info

! Command line arguments

character(rads_naml) :: filename
integer(fourbyteint) :: cyc, pass

! Formats

550 format (a)
551 format (a,' ... ',$)

! Initialise

call synopsis ('--head')
call rads_init (S)

! Load the surface_type grid

call getenv ('ALTIM', filename)
filename = trim(filename) // '/data/landmask.nc'
write (*,551) 'Loading mask '//trim(filename)
if (grid_load (filename, info) /= 0) call rads_exit ('Error loading landmask.')
write (*,550) 'done'

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

! Free the allocated grid

call grid_free(info)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Add surface type flags to RADS data', flag=flag)) return
call synopsis_devel ('')
write (*,1310)
1310  format (/ &
'This program changes BOTH the engineering flags and the surface_type variable')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: lon(n), lat(n), flags(n), surface_type(n)
integer(fourbyteint) :: flag, surf, i

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename)

! Get lat, lon, flags, surface_type

call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'flags', flags, .true.)
call rads_get_var (S, P, 'surface_type', surface_type, .true.)

! Process data records
!
! Bits and surface_type to be set when landmask contains any of the following values:
!
! landmask   | bit 4 | bit 5 | surface_type
! -----------------------------------------
! 0 = ocean  |   0   |   0   |   0
! 1 = land   |   1   |   1   |   2
! 2 = lake   |   0   |   1   |   3
! 3 = island |   1   |   1   |   2
! 4 = pond   |   0   |   1   |   3
! -----------------------------------------
!
! However, never undo surface_type = 4 (continental ice)

do i = 1,n

! Determine flags

	flag = nint(flags(i))
	surf = nint(surface_type(i))

	select case (nint(grid_query (info, lon(i), lat(i))))
	case (0) ! ocean
		flag = ibclr (flag, 4)
		flag = ibclr (flag, 5)
		surf = 0
	case (1, 3) ! land, island
		flag = ibset (flag, 4)
		flag = ibset (flag, 5)
		surf = 2
	case default ! lake, pond
		flag = ibclr (flag, 4)
		flag = ibset (flag, 5)
		surf = 3
	end select

	flags(i) = flag
	if (surface_type(i) /= 4) surface_type(i) = surf
enddo

! Store all data fields.

call rads_put_history (S, P)
call rads_def_var (S, P, 'flags')
call rads_def_var (S, P, 'surface_type')
call rads_put_var (S, P, 'flags', flags)
call rads_put_var (S, P, 'surface_type', surface_type)

write (*,552) n
end subroutine process_pass

end program rads_add_surface
