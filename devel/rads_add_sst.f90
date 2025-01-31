!-----------------------------------------------------------------------
! Copyright (c) 2011-2025  Remko Scharroo
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

!*rads_add_sst -- Add SST temperature and ice concentration to RADS data
!+
! This program adds the Reynolds SST and sea ice concentration
! to the RADS data. These fields are only provided as reference
! or can be used for data screening or validation purposed.
!
! The Reynolds SST grids and sea ice concentration grids are produced
! weekly and are converted first to NetCDF and stored in
! ${ALTIM}/data/sst. Note that the grids before 1990 were centered
! on Sundays, currently on Mondays.
!
! The SST grids have 1x1 degree resolution and contain no data gaps
! (invalid values or NaNs) and are interpolated linearly.
! Output units and resolution are the same as the input (0.01 degC).
!
! The sea ice concentation grids (1x1 deg) contain NaNs (over land);
! interpolation is linear with a minimum weight of 0.25.
! Output units and resolution are the same as the input (%).
!
! usage: rads_add_sst [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_sst

use rads
use rads_misc
use rads_devel
use rads_netcdf
use rads_time
use netcdf
use grib_api

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass
logical :: lice=.false., lsst=.false., lmean=.false., update=.false., new=.false.

! Data variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: t1old=0, t2old=0, j
integer(fourbyteint), parameter :: nx=1440, ny=720, secday=86400, sec2000=473299200
real(eightbytereal), parameter :: x0=0.125d0, y0=-89.875d0, dx=0.25d0, dy=0.25d0, dz=1d-2
integer(twobyteint) :: grids(0:nx+1,ny,2,3)

! Initialise

call synopsis ('--head')
call rads_set_options ('ismun ice sst mean all update new')
call rads_init (S)

! Get template for path name

call parseenv ('${ALTIM}/data/sst/oisst_v2.1/%Y%m/oisst-avhrr-v02r01.%Y%m%d.nc', path)

! Check all options
do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('i', 'ice')
		lice = .true.
	case ('s', 'sst')
		lsst = .true.
	case ('m', 'mean')
		lmean = .true.
	case ('all')
		lsst = .true.
		lmean = .true.
		lice = .true.
	case ('u', 'update')
		update = .true.
	case ('n', 'new')
		new = .true.
	end select
enddo
if (.not.lsst) update = .false.

! Process all data files

do cyc = S%cycles(1), S%cycles(2), S%cycles(3)
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cyc, pass, .true.)
		if (P%ndata > 0) call process_pass (P%ndata)
		call rads_close_pass (S, P)
	enddo
enddo

contains

!-----------------------------------------------------------------------
! Print synopsis
!-----------------------------------------------------------------------

subroutine synopsis (flag)
character(len=*), optional :: flag
if (rads_version ('Add SST temperature and ice concentration to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -i, --ice                 Add sea ice concentration' / &
'  -s, --sst                 Add sea surface temperature' / &
'  -m, --mean                Add local mean sea surface temperature' / &
'  --all                     All of the above' / &
'  -u, --update              Update files only when there are changes in SST (requires -s)' / &
'  -n, --new                 Only add sea ice when not yet existing (requires -u)')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, ix, iy, ncid
real(eightbytereal) :: time(n), lat(n), lon(n), ice(n), sst(n), ano(n), tmp(n)
integer(fourbyteint) :: t1, t2
real(eightbytereal) :: wx, wy, wt, f(2,2,2)
real(eightbytereal), parameter :: dz=1d-2
logical :: do_ice

call log_pass (P)

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)

! Process data records

do i = 1,n

! What was the start of the day?

	t1 = floor(time(i)/secday) * secday
	t2 = t1 + secday

! Load new grids when needed

	if (t1 == t2old) then ! Copy grids(t2) to grids(t1)
		grids(:,:,1,:) = grids(:,:,2,:)
		t1old = t1
	else if (t1 /= t1old) then ! Load new grids(t1)
		if (get_grid(t1,grids(:,:,1,:)) /= 0) then
			call log_string ('Warning: no SST for current time')
			call log_records (0)
			stop
		endif
		t1old = t1
	endif
	if (t2 /= t2old) then ! Load new grids(t2)
		if (get_grid(t2,grids(:,:,2,:)) /= 0) then
			call log_string ('Warning: no SST for current time')
			call log_records (0)
			stop
		endif
		t2old = t2
	endif

! Set weights for bi-linear interpolation in space

	if (lon(i) < -360d0 .or. lat(i) < -360d0) then
		sst(i) = nan
		ano(i) = nan
		ice(i) = nan
		cycle
	endif
	if (lon(i) < 0d0) lon(i) = lon(i) + 360d0
	wx = (lon(i)-x0)/dx + 1
	wy = (lat(i)-y0)/dy + 1
	ix = floor(wx)
	iy = floor(wy)
	wx = wx - ix
	wy = wy - iy

	f(1,:,:) = (1d0-wx)
	f(2,:,:) = wx
	f(:,1,:) = f(:,1,:) * (1d0-wy)
	f(:,2,:) = f(:,2,:) * wy

! Set weights for linear interpolation in time

	wt = (time(i) - t1)/(t2 - t1)
	f(:,:,1) = f(:,:,1) * (1d0-wt)
	f(:,:,2) = f(:,:,2) * wt

! Interpolate the global grids

	sst(i) = mat_product(grids(ix:ix+1,iy:iy+1,:,1),f)
	ano(i) = mat_product(grids(ix:ix+1,iy:iy+1,:,2),f)
	ice(i) = mat_product(grids(ix:ix+1,iy:iy+1,:,3),f)
enddo

! If requested, check for changes in sst first

if (update) then
	i = rads_verbose; rads_verbose = -1 ! Temporarily suspend warning
	call rads_get_var (S, P, 'sst', tmp, .true.)
	rads_verbose = i
	do i = 1,n
		if (isnan_(tmp(i)) .and. isnan_(sst(i))) cycle
		if (isnan_(tmp(i))) exit
		if (nint(tmp(i)/dz) /= nint(sst(i))) exit
	enddo
	if (i > n) then	! No changes
		call log_records (0)
		return
	endif
endif

! If "new" option is used, write seaice only when not yet existing

ncid = P%fileinfo(1)%ncid
do_ice = lice .and. .not.(new .and. nff(nf90_inq_varid(ncid,'seaice_conc',i)))

! Store all data fields

call rads_put_history (S, P)

if (do_ice) call rads_def_var (S, P, 'seaice_conc')
if (lsst)   call rads_def_var (S, P, 'sst')
if (lmean)  call rads_def_var (S, P, 'sst_mean')

if (do_ice) call rads_put_var (S, P, 'seaice_conc', ice)
if (lsst)   call rads_put_var (S, P, 'sst', sst*dz)
if (lmean)  call rads_put_var (S, P, 'sst_mean', (sst-ano)*dz)

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Get the SST and ice concentration grids for time t from NetCDF files
!-----------------------------------------------------------------------

integer function get_grid (t, grid)
integer(fourbyteint), intent(in) :: t
integer(twobyteint), intent(out) :: grid(0:nx+1,ny,3)
integer(fourbyteint) ::	ncid, varid, i, l, strf1985
character(len=rads_cmdl) :: fn

600 format ('(',a,')')

! Determine file name

l = strf1985 (fn, trim(path), t)
i = index(fn,'/',.true.)

! Open input file

if (nff(nf90_open(fn(:l), nf90_nowrite, ncid))) then
	write (*,600,advance='no') fn(i+1:l)
else
	! If not available, try preliminary
	write (*,600,advance='no') fn(i+1:l-3)//'_preliminary.nc'
	if (nft(nf90_open(fn(:l-3)//'_preliminary.nc', nf90_nowrite, ncid))) then
		get_grid = 1
		return
	endif
endif

! Get SST, SST anomaly and ice concentration grid

call nfs(nf90_inq_varid(ncid, 'sst', varid))
call nfs(nf90_get_var(ncid, varid, grid(1:nx,:,1)))
call nfs(nf90_inq_varid(ncid, 'anom', varid))
call nfs(nf90_get_var(ncid, varid, grid(1:nx,:,2)))
call nfs(nf90_inq_varid(ncid, 'ice', varid))
call nfs(nf90_get_var(ncid, varid, grid(1:nx,:,3)))

! Copy left and right boundaries

grid(0,:,:) = grid(nx,:,:)
grid(nx+1,:,:) = grid(1,:,:)

! Close file

call nfs(nf90_close(ncid))

get_grid = 0
end function get_grid

!-----------------------------------------------------------------------
! Weighted product with value checking on array a
!-----------------------------------------------------------------------

function mat_product (a, b)
real(eightbytereal) :: mat_product
integer(twobyteint), intent(in) :: a(8)
real(eightbytereal), intent(in) :: b(8)
real(eightbytereal) :: w
integer(fourbyteint) :: i
mat_product = 0d0
w = 0d0
do i = 1,8
	if (a(i) /= -999) then
		w = w + b(i)
		mat_product = mat_product + a(i)*b(i)
	endif
enddo
if (w < 0.5d0) then
	mat_product = 0d0
	w = 0d0
endif
mat_product = mat_product / w
end function mat_product

end program rads_add_sst
