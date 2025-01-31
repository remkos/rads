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

!*rads_add_era5 -- Add ERA5 meteo models to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! with values computed from ECMWF ERA5 meteorological models.
! The models provide sea level pressure, columnal water vapour content,
! surface temperature, 10-metre wind speed, and wave height.
!
! Input grids are found in the directory ${ALTIM}/data/era5.
!
! Interpolation is performed in hourly 0.5x0.5 degree grids;
! bi-linear in space, linear in time; contained in hourly files, with the
! file format:
! YYYY/MM/DD/ATIA_ANA_xx_xxx_<T-1H>_<T+1H>_<T>_ECMWF_ERA5.grib
! where:
! <T> is the epoch of the model in format YYYYMMDDHH0000Z
! <T-1H> and <T+1H> are the epoch minus 1 hour and plus 1 hour in the same format
!
! usage: rads_add_era5 [data-selectors] [options]
!
! References:
!
! Hopfield, H. S. (1969), Two-quartic tropospheric refractivity profile for
! correcting satellite data, J. Geophys. Res., 74(18), 4487-4499,
! 10.1029/JC074i018p04487.
!
! Saastamoinen, J. (1972), Atmospheric corrections for the troposphere and
! stratosphere in radio ranging of satellites, in The Use of Artificial
! Satellites for Geodesy, Geophys. Monogr. Ser., vol. 15, edited by
! S. W. Hendriksen, A. Mancini, and B. H. Chovitz, pp. 247-251,
! American Geophysical Union, Washington, D.C.
!
! Bevis, M., S. Businger, S. Chriswell, T. A. Herring, R. A. Anthes, C. Rocken,
! and R. H. Ware (1994), GPS meteorology: Mapping zenith wet delays onto
! precipitable water, J. Applied Meteorology, 33(3), 379-386,
! 10.1175/1520-0450(1994)033<0379:GMMZWD>2.0.CO;2.
!
! Mendes, V. B., G. Prates, L. Santos, and R. B. Langley (2000), An evaluation
! of the accuracy of models for the determination of the weighted mean temperature
! of the atmosphere, in Proc. of the 2000 National Technical Meeting of The
! Institute of Navigation, Anaheim, CA, January 2000, pp. 433-438.
!
! Petit, G., and B. Luzum (Eds.) (2010), IERS Conventions (2010),
! IERS Technical Note 36, Verlag des Bundesamts für Kartographie und Geodäsie.
!
! ECMWF (2022), ERA5, https://confluence.ecmwf.int/display/CKB/ERA5
!-----------------------------------------------------------------------
program rads_add_era5

use rads
use rads_misc
use rads_time
use rads_netcdf
use rads_devel
use netcdf
use meteo_subs
use grib_api

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: j, cyc, pass

! Data elements

character(rads_cmdl) :: pathnm
real(eightbytereal), parameter :: rad2=2d0*atan(1d0)/45d0
logical :: new = .false., test = .false., error

! Model data

integer, parameter :: mvar = 7, mext = mvar + 3
type :: var_
	logical requested, available
	character(len=80) :: shortName, name
end type
integer, parameter :: j_2t = 1, j_msl = 2, j_tcwv = 3, j_ib = mvar + 1, j_dry = mvar + 2, j_wet = mvar + 3

type(var_) :: var(mext)

integer(fourbyteint), parameter :: nx = 720, ny = 361, np = nx*ny
real(eightbytereal), parameter :: x0 = 0d0, x1 = 359.5d0, dx = 0.5d0
real(eightbytereal), parameter :: y0 = 90d0, y1 = -90d0, dy = -0.5d0
type :: model_
	integer(fourbyteint) :: hour
	real(eightbytereal) :: var(np,mvar)
end type
type(model_) :: m0, m1
integer(fourbyteint) :: fileid = -1, gribid = -1

! Model parameters, see Bevis et al [1994]

real(eightbytereal), parameter :: Rw = 461.5d3	! Gas constant water vapour per volume (Pa/K)
real(eightbytereal), parameter :: k2p = 0.221d0, k3 = 3739d0	! Refractivity constants (K/Pa, K^2/Pa)

! Initialise

var(j_2t)   = var_ (.false., .false., '2t',   '')
var(j_msl)  = var_ (.false., .false., 'msl',  '')
var(j_tcwv) = var_ (.false., .false., 'tcwv', '')
var(4)      = var_ (.false., .false., '10u',  'wind_speed_era5_u')
var(5)      = var_ (.false., .false., '10v',  'wind_speed_era5_v')
var(6)		= var_ (.false., .false., 'ci',   'seaice_conc_era5')
var(7)      = var_ (.false., .false., 'swh',  'swh_era5')
var(j_ib)   = var_ (.false., .false., '',     'inv_bar_static_era5')
var(j_dry)  = var_ (.false., .false., '',     'dry_tropo_era5')
var(j_wet)  = var_ (.false., .false., '',     'wet_tropo_era5')

m0%hour = -99999
m1%hour = -99999

call synopsis ('--head')
call rads_set_options ('n new all test')
call rads_init (S)

! Get template for path name

call parseenv ('${ALTIM}/data/era5/', pathnm)

! Which corrections are to be provided?

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('n', 'new')
		new = .true.
	case ('test')
		test = .true.
	case ('all')
		where (var(:)%name /= '') var(:)%requested = .true.
	end select
enddo

! Check requested output variables

do j = 1, S%nsel
	where (var(:)%name == S%sel(j)%name) var(:)%requested = .true.
enddo

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
if (rads_version ('Add ERA5 meteo models to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  --all                     ' &
'Same as -Vdry_tropo_era5,wet_tropo_era5,inv_bar_static_era5,swh_era5,sea_ice_conc_era5,', &
	'wind_speed_era5_u,wind_speed_era5_v' / &
'  -n, --new                 Only add variables when not yet existing' / &
'  --test                    Test/debug output')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, ncid
real(eightbytereal) :: time(n), lat(n), lon(n), h(n), surface_type(n), z(n, mext), &
	f0, f1, g1, g2, slp0, tmp
real(eightbytereal) :: w(4), x, y, wx, wy
integer :: ix, iy, idx(4), hour1, hour2, max_hour

! Formats

call log_pass (P)

! If 'new' option is used, write only when fields are not yet available

ncid = P%fileinfo(1)%ncid
if (new .and. nff(nf90_inq_varid(ncid,'dry_tropo_era5',i)) .and. &
	nff(nf90_inq_varid(ncid,'wet_tropo_era5',i))) then
	call log_records (0)
	return
endif

! Get time and position

call rads_get_var (S, P, 'time', time, .true.)
call rads_get_var (S, P, 'lat', lat, .true.)
call rads_get_var (S, P, 'lon', lon, .true.)
call rads_get_var (S, P, 'topo', h, .true.)
call rads_get_var (S, P, 'surface_type', surface_type, .true.)

! Correct DNSC08 or DTU10 topography of the Caspian Sea (-R46.5/54.1/36.5/47.1)
! It has bathymetry in both models in stead of lake topography (about -27 m)

do i = 1,n
	if (lon(i) > 46.5d0 .and. lon(i) < 54.1d0 .and. &
		lat(i) > 36.5d0 .and. lat(i) < 47.1d0 .and. &
		surface_type(i) < 2.5d0) h(i) = -27d0
enddo

! Get global pressure

call globpres(4,P%equator_time,slp0)

! Process data records

do i = 1,n
	x = time(i)/3600d0
	hour1 = floor(x)	! Counter of hourly periods
	hour2 = hour1 + 1
	max_hour = max(m0%hour, m1%hour)

! Load grid for first hour if both loaded grids are outdated

	if (hour1 > max_hour) error = get_gribs (hour1,m0)

! Load grid for second hour if both loaded grids are outdated

	if (hour2 > max_hour .and. .not.error) then
		if (hour1 == m0%hour) then
			error = get_gribs (hour2,m1)
		else
			error = get_gribs (hour2,m0)
		endif
	endif

! If error, escape

	if (error) then
		write (*,'(a)') 'Model switched off.'
		var(:)%requested = .false.
		exit
	endif

! Linearly interpolation in time, bi-linear interpolation in space

	f1 = abs(x - m0%hour)
	f0 = 1d0 - f1

! Find nearest grid points

	x = (lon(i)-x0)/dx
	ix = floor(x)
	wx = x - ix

	y = (lat(i)-y0)/dy
	iy = floor(y)
	wy = y - iy

	idx(1) = iy * nx + modulo(ix, nx) + 1
	idx(2) = iy * nx + modulo(ix+1, nx) + 1
	idx(3:4) = idx(1:2) + nx

	w(1) = (1d0-wx) * (1d0-wy)
	w(2) =      wx  * (1d0-wy)
	w(3) = (1d0-wx) *      wy
	w(4) =      wx  *      wy

! Interpolate all input variables

	do j=1,mvar
		if (var(j)%available) z(i,j) = f0 * dot_product (w,m0%var(idx,j)) + f1 * dot_product (w,m1%var(idx,j))
	enddo

! Correct sea level pressure (mbar) grids for altitude over land and lakes using Hopfield [1969]
! Convert pressure (mbar) to simple IB (m) if requested and over ocean
! Convert sea level pressure (mbar) to dry correction (m) using Saastamoinen models as referenced
! in IERS Conventions Chap 9

	if (surface_type(i) > 1.5d0) then ! Land and lakes
		g1 = -22768d-7 / (1d0 - 266d-5*cos(lat(i)*rad2) - 28d-8*h(i)) * pp_hop(h(i),lat(i),z(i,j_2t))
		g2 = 0d0	! Set IB to zero
	else    ! Ocean (altitude is zero)
		g1 = -22768d-7 / (1d0 - 266d-5*cos(lat(i)*rad2))
		g2 = -9.948d-3
	endif

! Derive dry tropospheric correction and IB

	if (var(j_dry)%requested .or. var(j_ib)%requested) then
		z(i,j_msl) = z(i,j_msl) * 1d-2 ! Convert from Pa to mbar

		! Convert sea level pressure to dry tropo correction after Saastamoinen [1972]
		z(i,j_dry) = z(i,j_msl) * g1

		! Convert sea level pressure to static inverse barometer
		z(i,j_ib) = (z(i,j_msl) - slp0) * g2
	endif

! Interpolate integrated water vapour and surface temperature in space and time

	if (var(j_wet)%requested) then
		! Convert surface temperature to mean temperature after Mendes et al. [2000]
		tmp = 50.4d0 + 0.789d0 * z(i,j_2t)
		! Convert integrated water vapour and mean temp to wet tropo correction
		! Also take into account conversion of wvc from kg/m^3 (= mm) to m.
		z(i,j_wet) = -1d-9 * Rw * (k3 / tmp + k2p) * z(i,j_tcwv)
	endif

	if (test) write (*,'(a26,2f6.3,2f12.6,2i4,6f7.3,4f9.4)') strf1985f(time(i),'T'), f0, f1, &
		lon(i), lat(i), ix+1, iy+1, wx, wy, w, z(i,j_tcwv), z(i,j_dry), z(i,j_wet), z(i,j_ib)
	if (test) write (*,'(8f9.4)') m0%var(idx,j_tcwv), m1%var(idx,j_tcwv)

enddo

! If no more fields are determined, abort.

if (.not.any(var(:)%requested)) then
	call log_records (0)
	stop
endif

! Store all data fields

call rads_put_history (S, P)

do j = 1,mext
	if (var(j)%requested) call rads_def_var (S, P, var(j)%name)
enddo

do j = 1,mext
	if (var(j)%requested) call rads_put_var (S, P, var(j)%name, z(:,j))
enddo

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! get_gribs -- Load necessary ECMWF meteo grids
!-----------------------------------------------------------------------
function get_gribs (hour, model)
integer(fourbyteint), intent(in) :: hour
type(model_), intent(inout) :: model
logical :: get_gribs
!
! Input are hourly files with all required fields
!--
! <hour> specifies the number of hours since 1 Jan 1985.
! Data is stored in a buffer <model>
!-----------------------------------------------------------------------
character(len=rads_cmdl) :: filenm
character(len=8) :: shortName
integer(fourbyteint) :: i,j,l,status,strf1985,dataTime,dataDate
real(eightbytereal) :: sec85

600 format (a)
1300 format (a,': ',a)

get_gribs = .true.

! Load new file and load common variables

if (test) write (*,*) 'hour = ',hour

l = strf1985(filenm, '%Y/%m/%d/ATIA_ANA_xx_xxx_00000000000000Z_00000000000000Z_%Y%m%d%H0000Z_ECMWF_ERA5.grib', &
	hour*3600)
l = strf1985(filenm(28:37), '%Y%m%d%H', (hour-1)*3600)
l = strf1985(filenm(44:53), '%Y%m%d%H', (hour+1)*3600)
l = len_trim(filenm)

write (*,600,advance='no') '(' // filenm(:l)

call grib_open_file(fileid,trim(pathnm) // filenm,'r',status)
if (status /= grib_success) then
	write (*,1300) 'Error opening file',trim(pathnm) // filenm(:l)
	return
endif

call grib_new_from_file(fileid,gribid,status)
if (status /= grib_success) then
	write (*,1300) 'Error loading first grib in',filenm(:l)
	return
endif

! Check the size
call grib_get(gribid,'numberOfPoints',i)
if (i /= np) then
	write (*,1300) 'numberOfPoints is not correct',filenm(:l)
	return
endif

! Check if we have the right dataDate and dataTime

call grib_get(gribid,'dataDate',dataDate)
call grib_get(gribid,'dataTime',dataTime)
dataTime = nint(sec85(2,dble(dataDate))/3600d0) + dataTime/100
if (dataTime /= hour) then
	write (*,1300) 'Incorrect date of',filenm(:l)
	return
endif

! Store the time

model%hour = hour

! Now read gribs
do
	call grib_get(gribid,'shortName',shortName)
	do j = 1,mvar
		if (var(j)%shortName == shortName) then
			write (*,600,advance='no') ' ' // trim(shortname)
			call grib_get(gribid,'values',model%var(:,j))
			var(j)%available = .true.
		endif
	enddo
	call grib_release(gribid)
	call grib_new_from_file(fileid,gribid,status)
	if (status /= grib_success) exit
enddo

! Close file
call grib_release(gribid)
call grib_close_file(fileid)
write (*,600,advance='no') ')'
if (test) write (*,600)

! Successful completion
get_gribs = .false.
end function get_gribs

end program rads_add_era5
