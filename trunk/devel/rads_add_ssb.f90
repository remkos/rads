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

!*rads_add_ssb -- Add (hybrid) SSB model to RADS data
!
! This program adds (or overwrites) the hybrid SSB model to RADS data.
! This is done by interpolating a grid in sigma0-SWH space, of which the
! filename is given in the <parameter> field in the rads.xml file.
!
! High and low SWH and sigma0 values are clipped to the extremes of the
! grid first.
!
! usage: rads_add_ssb [data-selectors]
!-----------------------------------------------------------------------
program rads_add_ssb

use rads
use rads_grid
use rads_misc
use rads_devel
use meteo_subs

! Data variables

type(rads_sat) :: S
type(rads_pass) :: P

! Other local variables

character(len=rads_cmdl) :: path
character(len=rads_varl) :: ssb_model = 'ssb_hyb'
integer(fourbyteint) :: i, cyc, pass
type(grid) :: info
type(rads_var), pointer :: var
logical :: lssb = .false., lwind = .false., saral

! Scan command line for options

call synopsis ('--head')
call rads_set_options ('s::w ssb:: wind all')
call rads_init (S)
saral = (S%sat == 'sa')
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('s', 'ssb')
		lssb = .true.
		if (rads_opt(i)%arg /= '') ssb_model = rads_opt(i)%arg(:rads_varl)
	case ('w', 'wind')
		lwind = .true.
	case ('all')
		lssb = .true.
		lwind = .true.
	end select
enddo

! Load SSB model

if (lssb) then
	var => rads_varptr (S, ssb_model)
	call parseenv ('${ALTIM}/data/models/' // var%info%parameters, path)
	if (grid_load(path,info) /= 0) call rads_exit ('Error loading '//trim(path))
endif

! Run process for all files

do cyc = S%cycles(1),S%cycles(2),S%cycles(3)
	do pass = S%passes(1),S%passes(2),S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

if (lssb) call grid_free (info)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('$Revision$', 'Add SSB model to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310 format (/ &
'Additional [processing_options] are:' / &
'  -s, --ssb[=MODEL]         Add/replace hybrid SSB model (default: ssb_hyb)' / &
'  -w, --wind                Compute altimeter wind speed from ECMWF model' / &
'  --all                     All of the above')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: x, y, xval(n), yval(n), ssb(n), wind(n)
integer :: i

! Formats

551 format (a,' ...',$)
552 format (i5,' records changed')

write (*,551) trim(P%filename(len_trim(S%dataroot)+2:))

! Compute wind speed after loading sigma0

if (lwind) then
	call rads_get_var (S, P, 'sig0', xval, .true.)
	wind = wind_ecmwf (xval, saral)
endif

! To compute SSB:
! load wind or sigma0 depending on the x-coordinate of the SSB grid, and load swh always

if (lssb) then
	if (index(info%xname, 'backscatter') > 0 .or. index(info%xname, 'sigma0') > 0) then
		if (.not.lwind) call rads_get_var (S, P, 'sig0', xval, .true.)
	else if (lwind) then
		xval = wind
	else
		call rads_get_var (S, P, 'wind_speed_alt', xval, .true.)
	endif
	call rads_get_var (S, P, 'swh', yval, .true.)
	do i = 1,n
		x = xval(i)
		if (x < info%xmin) x = info%xmin
		if (x > info%xmax) x = info%xmax
		y = yval(i)
		if (y < info%ymin) y = info%ymin
		if (y > info%ymax) y = info%ymax
		ssb(i) = grid_lininter (info, x, y)
	enddo
endif

! Write out all the data

call rads_put_history (S, P)

if (lwind) call rads_def_var (S, P, 'wind_speed_alt')
if (lssb)  call rads_def_var (S, P, ssb_model)

if (lwind) call rads_put_var (S, P, 'wind_speed_alt', wind)
if (lssb)  call rads_put_var (S, P, ssb_model, ssb)

write (*,552) n

end subroutine process_pass

end program rads_add_ssb
