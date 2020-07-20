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
use rads_misc
use rads_grid
use rads_devel

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P
type(grid) :: info
type(rads_var), pointer :: var
integer :: i

! Command line arguments

character(rads_cmdl) :: filename
integer(fourbyteint) :: cyc, pass
logical :: sentinel = .false.

! Initialise

call synopsis ('--head')
call rads_set_options ('s sentinel')
call rads_init (S)

! Check options

do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('s', 'sentinel')
		sentinel = .true.
	end select
enddo

! Load the surface type/class grid

if (sentinel) then
	var => rads_varptr (S, 'surface_class')
else
	var => rads_varptr (S, 'surface_type')
endif
call parseenv ('${ALTIM}/data/' // var%info%parameters, filename)

call log_string ('Loading mask '//filename)
if (grid_load (filename, info) /= 0) call rads_exit ('Error loading surface type grid')
call log_string ('done', .true.)

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
if (rads_version ('Add surface type flags to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -s, --sentinel            Use Sentinel-3 7-level mask (default is GMT landmask)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: lon(n), lat(n), flags(n), surface_type(n)
integer(fourbyteint) :: flag, surf, i

call log_pass (P)

! Get lat, lon, flags, surface_type

call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'flags', flags, .true.)
call rads_get_var (S, P, 'surface_type', surface_type, .true.)

! Process data records
!
! 1) When using the old-school GMT landmask
!
! Bits and surface_type to be set when the GMT landmask contains any of the following values:
!
! landmask   | bit 2 | bit 4 | bit 5 | surface_type
! -----------------------------------------
! 0 = ocean  |   0   |   0   |   0   |   0
! 1 = land   |   0   |   1   |   1   |   3
! 2 = lake   |   0   |   0   |   1   |   2
! 3 = island |   0   |   1   |   1   |   3
! 4 = pond   |   0   |   0   |   1   |   2
! cont. ice  |   1   |   1   |   1   |   4
! -----------------------------------------
!
! However, never undo surface_type = 4 or bit 2 = set (continental ice)
!
! 2) When using the Sentinel-3 surface type flags
!
! surf_class             | bit 2 | bit 4 | bit 5 | surface_type
! -------------------------------------------------------------
! 0 = ocean              |   0   |   0   |   0   |   0
! 1 = land               |   0   |   1   |   1   |   3
! 2 = continental water  |   0   |   0   |   1   |   2
! 3 = aquatic vegetation |   0   |   1   |   1   |   3
! 4 = cont. ice/snow     |   1   |   1   |   1   |   4
! 5 = floating ice       |   0   |   1   |   1   |   5
! 6 = salted basin       |   0   |   1   |   1   |   6
! 7 = undefined (= land) |   0   |   1   |   1   |   3
! -------------------------------------------------------------

do i = 1,n

! Determine flags

	! Get the values from array
	flag = nint(flags(i))
	surf = nint(surface_type(i))

	! When using 7-level Sentinel-3 mask
	if (sentinel) then
		surf = nint(grid_query (info, lon(i), lat(i)))

	! When using old-style GMT land mask
	else if (btest(flag,2) .or. surf == 4) then	! continental ice
		surf = 4
	else
		surf = nint(grid_query (info, lon(i), lat(i)))
		if (surf == 1) surf = 3 ! land
		if (surf == 4) surf = 2 ! pond
	endif

	! Clear or set bits
	! bit 2: continental ice
	if (surf /= 4) then
		flag = ibclr (flag, 2)
	else
		flag = ibset (flag, 2)
	endif
	! bit 4: water/dry
	if (surf == 0 .or. surf == 2) then
		flag = ibclr (flag, 4)
	else
		flag = ibset (flag, 4)
	endif
	! bit 5: ocean/land
	if (surf == 0) then
		flag = ibclr (flag, 4)
	else
		flag = ibset (flag, 4)
	endif

	! Store the values in array
	flags(i) = flag
	surface_type(i) = surf
enddo

! Store all data fields.

call rads_put_history (S, P)
call rads_def_var (S, P, 'flags')
call rads_def_var (S, P, var)
call rads_put_var (S, P, 'flags', flags)
call rads_put_var (S, P, var, surface_type)

call log_records (n)
end subroutine process_pass

end program rads_add_surface
