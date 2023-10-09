!-----------------------------------------------------------------------
! Copyright (c) 2011-2023  Remko Scharroo
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

!*rads_add_ww3_222 -- Add SWH from WaveWatch3 model to RADS data
!+
! This program adds SWH from the WaveWatch 3 model (version 2.22) to the
! RADS data base.
!
! The SWH values are stored in monthly 0.5x0.5 degree GRIB2 grids with
! "layers" at 3-hour intervals. These files can be found in
! ${ALTIM}/data/ww3_222.
!
! The grids contain invalid values, whose mask changes with time, because
! of ice cover. The grids are linearly interpolated in space and time
! and require a minimum weight of 0.5 for the value to be used.
!
! usage: rads_add_ww3_222 [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_ww3_222

use rads
use rads_misc
use rads_time
use rads_devel
use grib_api

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass
logical :: update = .false., lswh = .false.

! Data elements

character(rads_cmdl) :: path
integer(fourbyteint), parameter :: mt=249
real(eightbytereal), parameter :: dz=1d-2
integer(fourbyteint) :: j, nx, ny, nt, mjd, yy, mm, dd, yymm=0
real(eightbytereal) :: x0, x1, dx, y0, y1, dy, t0=0d0, dt
integer(twobyteint), allocatable :: grids(:,:,:)
integer(twobyteint) :: subgrid(2,2,2)

! Initialise

call synopsis ('--head')
call rads_set_options ('su swh all update')
call rads_init (S)

! Get template for path name

call parseenv ('${ALTIM}/data/ww3_222/multi_1.glo_30m.hs.%Y%m.grb2', path)

! Check all options
do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('s', 'swh')
		lswh = .true.
	case ('all')
		lswh = .true.
	case ('u', 'update')
		update = .true.
	end select
enddo

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

deallocate (grids,stat=j)

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add SWH from WaveWatch3 (v2.22) model to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -s, --swh                 Add WaveWatch3 (v2.22) SWH' / &
'  --all                     All of the above' / &
'  -u, --update              Update files only when there are changes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
real(eightbytereal) :: time(n), lat(n), lon(n), tmp(n), ww3(n), z, w, x, y, &
	f1, f2, f(2,2,2)
integer(fourbyteint) :: i, ix, iy, it

call log_pass (P)

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)

! Process data records

do i = 1,n

! Check if we rolled into a new month
! If so, load a whole month of gridded SWH

	mjd = floor(time(i)/86400d0)
	call mjd2ymd (mjd+46066, yy, mm, dd)
	if (yy*100+mm /= yymm) then
		if (get_ww3(mjd)) then
			call log_string ('Warning: no WW3 data for current time')
			call log_records (0)
			stop
		endif
		call ymd2mjd (yy, mm, 01, mjd)
		t0 = (mjd - 46066) * 86400d0
		yymm = yy*100+mm
	endif

! Skip data beyond latitude limits

	if (lat(i) < y0 .or. lat(i) > y1) then
		ww3(i) = nan
		cycle
	endif
	if (lon(i) < x0) lon(i) = lon(i) + 360d0

! Linearly interpolate in space and time

	f2 = (time(i)-t0)/dt + 1
	it = int(f2)

	if (it >= nt) then
		call log_string ('Warning: no WW3 data for current time')
		call log_records (0)
		stop
	endif

	f2 = f2 - it
	f1 = 1d0 - f2

	x = (lon(i)-x0)/dx + 1
	y = (lat(i)-y0)/dy + 1
	ix = int(x)
	iy = int(y)
	x = x - ix
	y = y - iy

	f(1,1,:) = (1d0-x)*(1d0-y)
	f(1,2,:) = (1d0-x)*(    y)
	f(2,1,:) = (    x)*(1d0-y)
	f(2,2,:) = (    x)*(    y)
	f(:,:,1) = f(:,:,1)*f1
	f(:,:,2) = f(:,:,2)*f2

	subgrid(1,:,:) = grids(ix,iy:iy+1,it:it+1)
	ix = mod(ix,nx)+1
	subgrid(2,:,:) = grids(ix,iy:iy+1,it:it+1)

	w = 0d0
	z = 0d0
	do it = 1,2
		do iy = 1,2
			do ix = 1,2
				if (subgrid(ix,iy,it) /= 32767) then
					w = w + f(ix,iy,it)
					z = z + f(ix,iy,it) * subgrid(ix,iy,it)
				endif
			enddo
		enddo
	enddo

	if (w > 0.5d0) then
		ww3(i) = z/w
	else
		ww3(i) = nan
	endif
enddo

! If requested, check for changes first

if (update) then
	i = rads_verbose; rads_verbose = -1 ! Temporarily suspend warning
	call rads_get_var (S, P, 'swh_ww3', tmp, .true.)
	rads_verbose = i
	do i = 1,n
		if (isnan_(tmp(i)) .and. isnan_(ww3(i))) cycle
		if (isnan_(tmp(i))) exit
		if (nint(tmp(i)/dz) /= nint(ww3(i))) exit
	enddo
	if (i > n) then	! No changes
		call log_records (0)
		return
	endif
endif

! Define and store data fields

call rads_put_history (S, P)
call rads_def_var (S, P, 'swh_ww3')
call rads_put_var (S, P, 'swh_ww3', ww3*dz)

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Get WaveWatch3 data
!-----------------------------------------------------------------------

function get_ww3 (mjd)
logical :: get_ww3
integer(fourbyteint), intent(in) :: mjd
integer(fourbyteint) ::	fileid, gribid, ix, iy, k, l, status, strf1985
real(eightbytereal), allocatable :: tmp(:)
character(len=rads_cmdl) :: fn
integer :: new(3), old(3)

600 format ('(',a,')')
1300 format (a,': ',a)

! Determine file name

get_ww3 = .true.
l = strf1985(fn, trim(path), mjd*86400)
write (*,600,advance='no') fn(l-29:l)

! Open input file

call grib_open_file(fileid,fn,'r',status)
if (status /= grib_success) then
	write (*,1300) 'Error opening file',fn(:l)
	return
endif
call grib_new_from_file(fileid,gribid,status)
if (status /= grib_success) then
	write (*,1300) 'Error loading grid',fn(:l)
	return
endif

! Get sizes

call grib_get(gribid,'Ni',nx,status)
call grib_get(gribid,'longitudeOfFirstGridPointInDegrees',x0,status)
call grib_get(gribid,'longitudeOfLastGridPointInDegrees',x1,status)
call grib_get(gribid,'iDirectionIncrementInDegrees',dx,status)

call grib_get(gribid,'Nj',ny,status)
call grib_get(gribid,'latitudeOfFirstGridPointInDegrees',y1,status)
call grib_get(gribid,'latitudeOfLastGridPointInDegrees',y0,status)
call grib_get(gribid,'jDirectionIncrementInDegrees',dy,status)

deallocate (grids,stat=k)
allocate (grids(nx,ny,mt),tmp(nx*ny))

! Get all grids

nt = 0
old = 0
do
	call grib_get(gribid,'dataDate',new(1))
	call grib_get(gribid,'dataTime',new(2))
	call grib_get(gribid,'forecastTime',new(3))
	if (all(new == old)) then
		write (*,1300) 'Warning: skipping duplicate message in',fn(:l)
		nt = nt -1
	endif
	old = new
	call grib_get(gribid,'values',tmp,status)
	call grib_release(gribid)
	if (status /= grib_success) exit
	nt = nt + 1
	k = 0
	do iy = ny,1,-1
		do ix = 1,nx
			k = k + 1
			grids(ix,iy,nt) = nint2(tmp(k)/dz)
		enddo
	enddo
	call grib_new_from_file(fileid,gribid,status)
	if (status /= grib_success) exit
enddo
deallocate(tmp)

! The time step varies between 3 and 6 hours. It depends on the number of
! time layers found. If more than 31*4 (124), then it must be 3 hours.

if (nt > 31*4) then
	dt = 3 * 3600d0
else
	dt = 6 * 3600d0
endif

! Close file

call grib_close_file(fileid)

get_ww3 = .false.
end function get_ww3

end program rads_add_ww3_222
