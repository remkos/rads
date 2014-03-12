!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2014  Remko Scharroo (Altimetrics LLC)
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

!*rads_fix_sa -- Patch RADS altimeter files of SARAL for various anomalies
!
! This program makes numerous patches to the SARAL RADS data processed
! by rads_gen_saral. These patches include:
!
! ssb:
! - Add/replace hybrid SSB to the RADS data
! wet:
! - Shift MWR wet prior to 2013-10-22: subtract 6.4 mm
! wind:
! - Update the wind speed using algorithm by Lillibrdge et al.
!
! usage: rads_fix_sa [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_fix_sa

use rads
use rads_misc
use rads_grid
use rads_devel
use meteo_subs
use netcdf
use rads_netcdf

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: i, cyc, pass
logical :: lssb = .false., lwet = .false., lwind = .false.
type(grid) :: issb_hyb

! Scan command line for options

call synopsis ('--head')
call rads_set_options (' sig0 ssb wet wind all')
call rads_init (S)
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('ssb')
		lssb = .true.
	case ('wet')
		lwet = .true.
	case ('wind')
		lwind = .true.
	case ('all')
		lssb = .true.
		lwet = .true.
		lwind = .true.
	end select
enddo

! Load SSB model

if (lssb) then
	call parseenv ('${ALTIM}/data/models/sa_ssb_hyb.nc?ssb_hyb', path)
	if (grid_load(path,issb_hyb) /= 0) call rads_exit ('Error loading '//trim(path))
endif

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	! Process passes
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

if (lssb) call grid_free (issb_hyb)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Patch SARAL data for several anomalies', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  --ssb                     Add/replace hybrid SSB model' / &
'  --wet                     Shift MWR wet prior to 2013-10-22 (pre-Patch2 data only)' / &
'  --wind                    Compute wind speed' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n),wet(n),swh(n),ssb(n),sig0(n),wind(n),x,y
real(eightbytereal), parameter :: time0 = 909058527d0 ! 2013-10-22 12:15:27

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename(len_trim(S%dataroot)+2:))

! Process data records

call rads_get_var (S, P, 'sig0_ka', sig0, .true.)

! Compute SSB

if (lssb) then
	call rads_get_var (S, P, 'swh_ka', swh, .true.)
	do i = 1,n
		x = sig0(i)
		if (x < issb_hyb%xmin) x = issb_hyb%xmin
		if (x > issb_hyb%xmax) x = issb_hyb%xmax
		y = swh(i)
		if (y < issb_hyb%ymin) y = issb_hyb%ymin
		if (y > issb_hyb%ymax) y = issb_hyb%ymax
		ssb(i) = grid_lininter (issb_hyb, x, y)
	enddo
endif

! Shift MWR wet tropo prior to 2013-10-22 12:15:27 (pre-patch2 data only)

if (P%start_time > time0 .or. index(P%original, '(V5') > 0) lwet = .false.
if (lwet) then
	call rads_get_var (S, P, 'time', time, .true.)
	call rads_get_var (S, P, 'wet_tropo_rad', wet, .true.)
	where (time < time0) wet = wet - 6.4d-3
endif

! Compute wind speed

if (lwind) wind = wind_ecmwf (sig0, .true.)

! If nothing changed, stop here

if (.not.(lssb .or. lwet .or. lwind)) then
	write (*,552) 0
	return
endif

! Write out all the data

call rads_put_history (S, P)
if (lssb) call rads_def_var (S, P, 'ssb_hyb')
if (lssb) call rads_put_var (S, P, 'ssb_hyb', ssb)
if (lwet) call rads_put_var (S, P, 'wet_tropo_rad', wet)
if (lwind) call rads_put_var (S, P, 'wind_speed_alt', wind)

write (*,552) n
end subroutine process_pass

end program rads_fix_sa
