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

!*rads_add_mog2d -- Add MOG2D dynamic atmospheric correction to RADS data
!+
! This program adds the MOG2D correction to RADS data. The MOG2D
! serves as a replacement to the conventional IB correction.
! The MOG2D model provides the high-frequency effect of wind and
! pressure loading on the ocean and is complemented by the
! low-frequency IB correction. The field number 1004 added to the
! data by this routine includes both the low- and high-frequency
! components.
!
! Input grids are found in the directory ${ALTIM}/data/dac
! The grids are bzip2 compressed NetCDF files. The required grids
! are unpacked and placed in /tmp before being used in this program.
!
! Interpolation is performed in 6-hourly grids of 0.25x0.25 degree
! spacing; bi-linear in space, linear in time. Note that the MOG2D
! grids have the order of the dimensions longitude and latitude
! flipped around from what is standard, i.e., (lon,lat) instead of
! (lat,lon).
!
! usage: rads_add_mog2d [data-selectors] [options]
!-----------------------------------------------------------------------
program rads_add_mog2d

use rads
use rads_misc
use rads_devel
use rads_netcdf
use netcdf

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint), parameter :: type_mog2d=0, type_era=1
integer(fourbyteint) :: cyc, pass, j, type = type_mog2d
logical :: update = .false.


! Data variables

character(rads_cmdl) :: path
character(rads_varl) :: varnm
integer(fourbyteint) :: hexold=-99999, first=1
integer(fourbyteint), parameter :: nx=1441, ny=721
real(eightbytereal), parameter :: x0=0d0, y0=-90d0, dx=0.25d0, dy=0.25d0
real(eightbytereal) :: z0, dz
integer(twobyteint) :: grids(nx,ny,2), tmp(ny,nx-1)

! Initialise

call synopsis ('--head')
call rads_set_options ('uet update era all')
call rads_init (S)

! Check all options
do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('u', 'update')
		update = .true.
	case ('e', 'era')
		type = type_era
	end select
enddo

! Get template for path name and set variable name

select case (type)
case (type_era)
	call parseenv ('${RADSROOT}/ext/slcci/Products/DAC_ERA_Interim_20161115/%Y/dac_era_interim_%Y%m%d_%H.nc', path)
	varnm = 'inv_bar_mog2d_era'
case default
	call parseenv ('${ALTIM}/data/dac/%Y/dac_dif_%Y%m%d_%H.nc', path)
	varnm = 'inv_bar_mog2d'
end select

! Link variable


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
if (rads_version ('Add MOG2D dynamic atmospheric correction to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  --all                     (Has no effect)'/ &
'  -e, --era                 Use the DAC forced by ERA Interim (1991-2015 only)'/ &
'  -u, --update              Update files only when there are changes')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, ix, iy, hex
real(eightbytereal) :: time(n), lat(n), lon(n), cor(n), tmp(n)
real(eightbytereal) :: f1, f2, f(2,2), z1, z2, x, y
logical :: err

call log_pass (P)

! Get time and location

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)

! Process data records

do i = 1,n

	f2 = time(i)/21600d0 ! Number of 6-hourly blocks since 1985
	hex = int(f2)

! Load new grids when entering new 6-hour period

	if (hex /= hexold) then
		if (hex == hexold+1) then
			! Replace first (oldest) grid with a new grid
			err = get_mog2d(hex+1,grids(:,:,first))
			first = 3-first	! Switch notion of first and second grid
		else
			! Replace both grids
			first = 1
			err = get_mog2d(hex,grids(:,:,1)) .or. get_mog2d(hex+1,grids(:,:,2))
		endif
		if (err) then
			call log_string ('Warning: No MOG2D field for current time')
			call log_records (0)
			stop
		endif
		hexold = hex
	endif

! Linearly interpolate in space and time

	if (lon(i) < 0d0) lon(i) = lon(i) + 360d0
	f2 = f2 - hex
	f1 = 1d0 - f2

	x = (lon(i)-x0)/dx+1
	y = (lat(i)-y0)/dy+1
	ix = int(x)
	iy = int(y)
	x = x - ix
	y = y - iy

	f(1,1) = (1d0-x)*(1d0-y)
	f(1,2) = (1d0-x)*(    y)
	f(2,1) = (    x)*(1d0-y)
	f(2,2) = (    x)*(    y)

	z1 = sum(grids(ix:ix+1,iy:iy+1,  first)*f)
	z2 = sum(grids(ix:ix+1,iy:iy+1,3-first)*f)

	cor(i) = (f1*z1 + f2*z2) * dz + z0
enddo

! If requested, check for changes first

if (update) then
	i = rads_verbose; rads_verbose = -1 ! Temporarily suspend warning
	call rads_get_var (S, P, varnm, tmp, .true.)
	rads_verbose = i
	do i = 1,n
		if (isnan_(tmp(i)) .and. isnan_(cor(i))) cycle
		if (isnan_(tmp(i))) exit
		if (nint(tmp(i)/dz) /= nint(cor(i)/dz)) exit
	enddo
	if (i > n) then	! No changes
		call log_records (0)
		return
	endif
endif

! Store all data fields

call rads_put_history (S, P)
call rads_def_var (S, P, varnm)
call rads_put_var (S, P, varnm, cor)

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! Get the MOG2D grid for "hex"
!-----------------------------------------------------------------------

function get_mog2d (hex, grid)
integer(fourbyteint) :: hex
logical :: get_mog2d
integer(twobyteint) :: grid(:,:)
character(len=rads_cmdl) :: filenm
integer(fourbyteint) ::	ncid,v_id,j,l,strf1985
real(eightbytereal) :: time

600 format ('(',a,')')
1300 format (a,': ',a)

get_mog2d = .true.

! Determine file name

l = strf1985(filenm, trim(path), hex*21600)

! Open input file

l = len_trim(filenm)
j = index(filenm, '/', .true.)
write (*,600,advance='no') filenm(j+1:l)
if (nft(nf90_open(filenm,nf90_nowrite,ncid))) then
	write (*,1300) 'Error opening file',filenm(:l)
	return
endif

! Check if NetCDF file contains variable name Grid_0001 or dac

if (nft(nf90_inq_varid(ncid,'dac',v_id))) then
	if (nft(nf90_inq_varid(ncid,'Grid_0001',v_id))) call fin('Error finding variable')
endif

! Get CNES JD and check against input

if (nft(nf90_get_att(ncid,v_id,'Date_CNES_JD',time))) call fin('Error reading time')
if (nint((time-12784)*24) /= hex*6) call fin('Hour does not match grid')

! Get scale factor, offset and missing value

if (nft(nf90_get_att(ncid,v_id,'add_offset',z0))) z0=0
if (nft(nf90_get_att(ncid,v_id,'scale_factor',dz))) dz=1

! Check if the grid is transposed

if (nft(nf90_inq_varid(ncid,'LatLon',j))) then
	! No LatLon variable: read normal grid
	if (nft(nf90_get_var(ncid,v_id,grid(1:nx-1,:)))) call fin('Error reading data grid')
else
	! If transposed, read temporary grid, then transpose
	if (nft(nf90_get_var(ncid,v_id,tmp))) call fin('Error reading data grid')
	grid(1:nx-1,:) = transpose(tmp)
endif

! Copy Greenwich meridian

grid(nx,:) = grid(1,:)

j = nf90_close(ncid)
get_mog2d = .false.
end function get_mog2d

end program rads_add_mog2d
