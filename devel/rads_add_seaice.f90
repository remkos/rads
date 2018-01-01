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

!*rads_add_seaice -- Add OSIAF sea ice concentration to RADS data
!+
! This program adds the Ocean and Sea Ice Satellite Application
! Facility (OSI SAF) sea ice concentration to the RADS data. This field
! is only provided as reference or can be used for data screening or
! validation purposed.
!
! The OSI SAF sea ice concentration grids are produced daily. In the
! PO.DAAC/NCEI version of GHRSST, the original OSI SAF grids have been
! interpolated from a polar stereographic projection to equally spaced,
! fixed grids. The netCDF files are stored in ${ALTIM}/data/OSAIF/[doy]
!
! The sea ice grids have 0.01x0.01 degree resolution. Both land and
! ice-free ocean locations contain NaNs; interpolation is linear with a
! minimum weight of 0.25. Output units and resolution are the same as
! the input (%).
!
! usage: rads_add_seaice [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_seaice

use rads
use rads_misc
use rads_devel
use rads_netcdf
use rads_time
use netcdf

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: cyc, pass
logical :: update=.false.

! Data variables

character(len=rads_cmdl) :: path
integer(fourbyteint) :: day1985old=-99999, first=1, j
integer(fourbyteint), parameter :: nx=36000, ny=18000
!real(eightbytereal) :: dz
real(eightbytereal), parameter :: x0=-180.00d0, y0=-89.99d0, dx=1d-2, dy=1d-2
integer(onebyteint) :: grids(nx,ny,2)

! Initialise

call synopsis ('--head')
call rads_set_options ('update')
call rads_init (S)

! Check all options
do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('u', 'update')
		update = .true.
	end select
enddo

! Get template for path name

call parseenv ('${ALTIM}/data/OSIAF/%Y/%j/%Y%m%d090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc',path)

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
if (rads_version ('Add  ice concentration to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -u, --update              Update files only when there are changes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, ix, iy, day1985
real(eightbytereal) :: time(n), lat(n), lon(n), ice(n), tmp(n)
integer(fourbyteint) :: t1, t2
real(eightbytereal) :: wx, wy, wt, f(2,2,2), f2
real(eightbytereal), parameter :: dz=1d-2
logical :: err

call log_pass (P)

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)

! Process data records

do i = 1,n

! Today and yesterday
   t1 = floor(time(i)/86400)
   t2 = t1 + 1

! Load new grids when needed

	f2 = time(i)/86400d0 ! Number of days since 1985
	day1985 = int(f2)

! Load new grids when entering new day

	if (day1985 /= day1985old) then
		if (day1985 == day1985old+1) then
			! Replace first (oldest) grid with a new grid
			err = get_osiaf(day1985+1,grids(:,:,first))
			first = 3-first	! Switch notion of first and second grid
		else
			! Replace both grids
			first = 1
			err = get_osiaf(day1985,grids(:,:,1)) .or. get_osiaf(day1985+1,grids(:,:,2))
		endif
		if (err) then
			call log_string ('Warning: No OSIAF field for current time')
			call log_records (0)
			stop
		endif
		day1985old = day1985
	endif

! Set weights for bi-linear interpolation in space

	if (lon(i) < -360d0 .or. lat(i) < -360d0) then
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

	wt = (time(i)/86400d0 - t1)/(t2 - t1)
	f(:,:,1) = f(:,:,1) * (1d0-wt)
	f(:,:,2) = f(:,:,2) * wt

! Interpolate ice concentration (has invalid values)

	ice(i) = mat_product(grids(ix:ix+1,iy:iy+1,:),f)
enddo

! If requested, check for changes in sea ice first

if (update) then
	i = rads_verbose; rads_verbose = -1 ! Temporarily suspend warning
	call rads_get_var (S, P, 'seaice_conc_osiaf', tmp, .true.)
	rads_verbose = i
	do i = 1,n
		if (isnan_(tmp(i)) .and. isnan_(ice(i))) cycle
		if (isnan_(tmp(i))) exit
		if (nint(tmp(i)/dz) /= nint(ice(i))) exit
	enddo
	if (i > n) then	! No changes
		call log_records (0)
		return
	endif
endif

! Store all data fields

call rads_put_history (S, P)

call rads_def_var (S, P, 'seaice_conc_osiaf')

call rads_put_var (S, P, 'seaice_conc_osiaf', ice)

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Weighted product with value checking on array a
!-----------------------------------------------------------------------

function mat_product (a, b)
real(eightbytereal) :: mat_product
integer(onebyteint), intent(in) :: a(8)
real(eightbytereal), intent(in) :: b(8)
real(eightbytereal) :: w
integer(fourbyteint) :: i
mat_product = 0d0
w = 0d0
do i = 1,8
	if (a(i) <= 100) then
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

!-----------------------------------------------------------------------
! Get the OSIAF sea ice concentration grid for "day1985"
!-----------------------------------------------------------------------

function get_osiaf (day1985, grid)
integer(fourbyteint) :: day1985
logical :: get_osiaf
integer(onebyteint) :: grid(:,:)
character(len=rads_cmdl) :: filenm
integer(fourbyteint) ::	ncid,v_id,t_id,j,l,strf1985
real(eightbytereal) :: time
real(eightbytereal) :: z0, dz

600 format ('(',a,')')
1300 format (a,': ',a)

get_osiaf = .true.

! Determine file name

l = strf1985(filenm, path, day1985*86400)

! Open input file

write (*,600,advance='no') trim(filenm)
if (nft(nf90_open(filenm,nf90_nowrite,ncid))) then
	write (*,1300) 'Error opening file',filenm(:l)
	return
endif

! Check if netCDF file contains variable name time

if (nft(nf90_inq_varid(ncid,'time',t_id))) call fin('Error finding variable')

! Get time in 1980 seconds and check against input

if (nft(nf90_get_var(ncid,t_id,time))) call fin('Error reading time')
if (nint((time-126230400)/86400) /= day1985) call fin('Day does not match grid')

! Check if netCDF file contains variable name sea_ice_fraction

if (nft(nf90_inq_varid(ncid,'sea_ice_fraction',v_id))) call fin('Error finding variable')

! Get scale factor, offset and missing value

if (nft(nf90_get_att(ncid,v_id,'add_offset',z0))) z0=0
if (nft(nf90_get_att(ncid,v_id,'scale_factor',dz))) dz=1

if (nft(nf90_get_var(ncid,v_id,grid(2:nx,:),count = (/ nx, ny-1, 1 /)))) call fin('Error reading data grid')

! Set missing value to 0

where (grid .eq. -128) grid = 0

! Copy Date Line meridian

grid(1,:) = grid(nx,:)

j = nf90_close(ncid)
get_osiaf = .false.
end function get_osiaf

end program rads_add_seaice
