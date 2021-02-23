!-----------------------------------------------------------------------
! Copyright (c) 2011-2021  Remko Scharroo
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

!*rads_add_era -- Add ERA interim meteo models to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! with values computed from ECMWF ERA interim meteorological models.
! The models provide sea level pressure, columnal water vapour content and
! surface temperature.
!
! Input grids are found in the directory ${ALTIM}/data/era-int-full.
!
! Interpolation is performed in 6-hourly reduced Gaussian grids (N128);
! bi-linear in space, linear in time; contained in monthly files.
!
! usage: rads_add_era [data-selectors] [options]
!
! Note: ECMWF stopped providing ERA interim model grids beyond August 2019,
! favouring ERA5 instead. A new programme was set up to process those
! model grids: rads_add_era5.
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
!-----------------------------------------------------------------------
program rads_add_era

use rads
use rads_misc
use rads_netcdf
use rads_devel
use netcdf
use meteo_subs
use tides
use grib_api

! Command line arguments

type(rads_sat) :: S
type(rads_pass) :: P
integer(fourbyteint) :: j, cyc, pass

! Data elements

character(rads_cmdl) :: path
integer(fourbyteint) :: hex, hexold=-99999
type(airtideinfo) :: airinfo
real(eightbytereal), parameter :: rad2=2d0*atan(1d0)/45d0
logical :: dry_on=.false., wet_on=.false., ib_on=.false., air_on=.false., new=.false., error

! Model data

integer(fourbyteint), parameter :: np_max=88838, nj_max=259
type :: model_
	integer(fourbytereal) :: nj, np, ip, pl(nj_max), ql(nj_max), idx(4)
	real(eightbytereal) :: glat(np_max), glon(np_max), w(4)
	real(eightbytereal) :: slp(np_max), wvc(np_max), tmp(np_max)
end type
type(model_) :: m1, m2

! Model parameters, see Bevis et al [1994]

real(eightbytereal), parameter :: Rw = 461.5d3	! Gas constant water vapour per volume (Pa/K)
real(eightbytereal), parameter :: k2p = 0.221d0, k3 = 3739d0	! Refractivity constants (K/Pa, K^2/Pa)

! Initialise

call synopsis ('--head')
call rads_set_options ('dwain dry wet air ib all new')
call rads_init (S)

! Get template for path name

call parseenv ('${ALTIM}/data/era-int-full/era-int.%Y%m.grb', path)

! Which corrections are to be provided?

do j = 1,rads_nopt
	select case (rads_opt(j)%opt)
	case ('d', 'dry')
		dry_on = .true.
	case ('w', 'wet')
		wet_on = .true.
	case ('a', 'air')
		air_on = .true.
	case ('i', 'ib')
		ib_on = .true.
	case ('n', 'new')
		new = .true.
	case ('all')
		dry_on = .true.
		wet_on = .true.
		air_on = .true.
		ib_on = .true.
	end select
enddo

! Init air tide

if (dry_on .or. ib_on) call airtideinit ('airtide', airinfo)

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
if (rads_version ('Add ERA interim meteo models to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -d, --dry                 Add ERA interim dry tropospheric correction' / &
'  -w, --wet                 Add ERA interim wet tropospheric correction' / &
'  -a, --air                 Add air tide' / &
'  -i, --ib                  Add static inverse barometer correction' / &
'  --all                     All of the above' / &
'  -n, --new                 Only add variables when not yet existing')
stop
end subroutine synopsis

!-----------------------------------------------------------------------
! Process a single pass
!-----------------------------------------------------------------------

subroutine process_pass (n)
integer(fourbyteint), intent(in) :: n
integer(fourbyteint) :: i, ncid
real(eightbytereal) :: time(n), lat(n), lon(n), h(n), surface_type(n), dry(n), wet(n), ib(n), air(n), &
	f1, f2, g1, g2, slp, dslp, slp0, wvc, tmp

! Formats

call log_pass (P)

! If "new" option is used, write only when fields are not yet available

ncid = P%fileinfo(1)%ncid
if (new .and. nff(nf90_inq_varid(ncid,'dry_tropo_era',i)) .and. &
	nff(nf90_inq_varid(ncid,'wet_tropo_era',i))) then
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
	f2 = time(i)/21600d0
	hex = floor(f2)	! Counter of 6-hourly periods

! Load new grids when entering new 6-hour period

	if (hex /= hexold) then
		if (hex == hexold + 1) then
			m1 = m2
			error = get_gribs (hex+1,m2)
		else
			error = get_gribs (hex,m1)
			if (.not.error) error = get_gribs (hex+1,m2)
		endif
		if (error) then
			write (*,'(a)') 'Model switched off.'
			dry_on = .false.
			ib_on = .false.
			wet_on = .false.
			exit
		endif
		hexold = hex
	endif

! Linearly interpolation in time, bi-linear interpolation in space

	if (lon(i) < 0d0) lon(i) = lon(i) + 360d0
	f2 = f2 - hex
	f1 = 1d0 - f2

! Find nearest grid points

	call find_nearest (lat(i), lon(i), m1)
	call find_nearest (lat(i), lon(i), m2)

! Interpolate surface temperature in space and time

	tmp = f1 * dot_product (m1%w,m1%tmp(m1%idx)) + f2 * dot_product (m2%w,m2%tmp(m2%idx))

! Correct sea level pressure grids for altitude over land and lakes using Hopfield [1969]
! Convert pressure (mbar) to simple IB (m) if requested and over ocean
! Convert sea level pressure (mbar) to dry correction (m) using Saastamoinen models as referenced
! in IERS Conventions Chap 9

	if (surface_type(i) > 1.5d0) then ! Land and lakes
		g1 = -22768d-7 / (1d0 - 266d-5*cos(lat(i)*rad2) - 28d-8*h(i)) * pp_hop(h(i),lat(i),tmp)
		g2 = 0d0	! Set IB to zero
	else    ! Ocean (altitude is zero)
		g1 = -22768d-7 / (1d0 - 266d-5*cos(lat(i)*rad2))
		g2 = -9.948d-3
	endif

! Interpolate sea level pressure in space and time and add airtide correction

	if (dry_on .or. ib_on .or. air_on) then
		slp = f1 * dot_product (m1%w,m1%slp(m1%idx)) + f2 * dot_product (m2%w,m2%slp(m2%idx))

		! Remove-and-restore the air tide
		dslp = airtide (airinfo, time(i), lat(i), lon(i)) &
			- f1 * airtide (airinfo, hex * 21600d0, lat(i), lon(i)) &
			- f2 * airtide (airinfo, (hex+1) * 21600d0, lat(i), lon(i))
		slp = slp * 1d-2 + dslp ! Convert Pa to hPa (mbar)

		! Convert sea level pressure to dry tropo correction after Saastamoinen [1972]
		dry(i) = slp * g1
		air(i) = dslp * g1

		! Convert sea level pressure to static inverse barometer
		ib(i) = (slp - slp0) * g2
	endif

! Interpolate integrated water vapour and surface temperature in space and time

	if (wet_on) then
		wvc = f1 * dot_product (m1%w,m1%wvc(m1%idx)) + f2 * dot_product (m2%w,m2%wvc(m2%idx))
		! Convert surface temperature to mean temperature after Mendes et al. [2000]
		tmp = 50.4d0 + 0.789d0 * tmp
		! Convert integrated water vapour and mean temp to wet tropo correction
		! Also take into account conversion of wvc from kg/m^3 (= mm) to m.
		wet(i) = -1d-9 * Rw * (k3 / tmp + k2p) * wvc
	endif
enddo

! If no more fields are determined, abort.

if (.not.(dry_on .or. ib_on .or. wet_on)) then
	call log_records (0)
	stop
endif

! Store all data fields

call rads_put_history (S, P)

if (dry_on) call rads_def_var (S, P, 'dry_tropo_era')
if (wet_on) call rads_def_var (S, P, 'wet_tropo_era')
if (ib_on ) call rads_def_var (S, P, 'inv_bar_static')
if (air_on) call rads_def_var (S, P, 'dry_tropo_airtide')

if (dry_on) call rads_put_var (S, P, 'dry_tropo_era', dry)
if (wet_on) call rads_put_var (S, P, 'wet_tropo_era', wet)
if (ib_on ) call rads_put_var (S, P, 'inv_bar_static', ib)
if (air_on) call rads_put_var (S, P, 'dry_tropo_airtide', air)

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! get_gribs -- Load necessary ECMWF meteo grids
!-----------------------------------------------------------------------
function get_gribs (hex, model)
integer(fourbyteint), intent(in) :: hex
type(model_), intent(inout) :: model
logical :: get_gribs
! Input are monthly files with all three required fields of the form:
! ${ALTIM}/data/era-int-full/era-int.201201.grb
!
! <hex> specifies the number of 6-hourly blocks since 1 Jan 1985.
! Data is stored in a buffer <model>
!-----------------------------------------------------------------------
character(len=rads_cmdl) :: filenm,old_filenm=''
character(len=8) :: shortName
integer(fourbyteint) :: fileid=-1,gribid=-1,i,l,status,strf1985,dataTime,dataDate
real(eightbytereal) :: sec85

600 format ('(',a,1x,i0,')')
1300 format (a,': ',a)

get_gribs = .true.

! Load new file and load common variables

l = strf1985(filenm, path, hex*21600)

if (filenm /= old_filenm) then
	if (gribid /= -1) call grib_release(gribid)
	if (fileid /= -1) call grib_close_file(fileid)
	write (*,600,advance='no') filenm(l-17:l),hex

	call grib_open_file(fileid,filenm,'r',status)
	if (status /= grib_success) then
		write (*,1300) 'Error opening file',filenm(:l)
		return
	endif

	call grib_new_from_file(fileid,gribid,status)
	if (status /= grib_success) then
		write (*,1300) 'Error loading first grib in',filenm(:l)
		return
	endif

	! Check sizes and get longitudes and latitudes
	call grib_get(gribid,'numberOfPoints',model%np)
	if (model%np > np_max) then
		write (*,1300) 'numberOfPoints is too large in',filenm(:l)
		return
	endif
	call grib_get(gribid,'Nj',model%nj)
	if (model%nj > nj_max) then
		write (*,1300) 'Nj is too large in',filenm(:l)
		return
	endif
	model%ip = 1
	call grib_get(gribid,'pl',model%pl)
	model%ql(1) = 1
	do i = 2,model%nj
		model%ql(i) = model%ql(i-1) + model%pl(i-1)
	enddo
	call grib_get_data(gribid,model%glat,model%glon,model%slp)
endif

! Loop until we find the right dataDate and dataTime
do
	call grib_get(gribid,'dataDate',dataDate)
	call grib_get(gribid,'dataTime',dataTime)
	dataTime = nint(sec85(2,dble(dataDate))/21600d0) + dataTime/600
	if (dataTime == hex) exit
	call grib_release(gribid)
	call grib_new_from_file(fileid,gribid,status)
	if (status /= grib_success) then
		write (*,1300) 'Reached end of',filenm(:l)
		return
	endif
enddo

! Now read gribs
call grib_get(gribid,'shortName',shortName)
if (shortName /= 'tcwv') then
	write (*,1300) 'Wrong field order ('//trim(shortName)//' != tcwv) in',filenm(:l)
	return
endif
call grib_get(gribid,'values',model%wvc)
call grib_release(gribid)
call grib_new_from_file(fileid,gribid,status)
call grib_get(gribid,'shortName',shortName)
if (shortName /= 'msl') then
	write (*,1300) 'Wrong field order ('//trim(shortName)//' != msl) in',filenm(:l)
	return
endif
call grib_get(gribid,'values',model%slp)
call grib_release(gribid)
call grib_new_from_file(fileid,gribid,status)
call grib_get(gribid,'shortName',shortName)
if (shortName /= '2t') then
	write (*,1300) 'Wrong field order ('//trim(shortName)//' != 2t) in',filenm(:l)
	return
endif
call grib_get(gribid,'values',model%tmp)

! Successful completion
get_gribs = .false.
end function get_gribs

!-----------------------------------------------------------------------
! find_nearest - Find the four grid points nearest to given point
!-----------------------------------------------------------------------

subroutine find_nearest (lat, lon, model)
real(eightbytereal), intent(in) :: lat, lon
type(model_), intent(inout) :: model
! Upon return model%idx and model%w are index and weight of closest points (NW, NE, SW, SE)
integer(fourbyteint) :: i,j,ip
real(eightbytereal) :: d

! Find highest latitude index for which glat >= lat; ip is north of lat, ip+1 is south of lat
ip = model%ip
do while (model%glat(model%ql(ip+1)) >= lat)
	ip = ip + 1
enddo
do while (model%glat(model%ql(ip  )) < lat)
	ip = ip - 1
enddo
model%ip = ip

! Compute longitude indices and weights at northern side (j=0) and southern side (j=1)
do j = 0,1
	d = lon / 360d0 * model%pl(ip+j)
	i = floor(d)
	d = d - i
	model%idx(2*j+1) = model%ql(ip+j) + modulo(i  ,model%pl(ip+j))
	model%idx(2*j+2) = model%ql(ip+j) + modulo(i+1,model%pl(ip+j))
	model%w(2*j+1) = 1d0 - d
	model%w(2*j+2) = d
enddo

! Multiply with weight in latitude direction
d = (lat - model%glat(model%ql(ip))) / (model%glat(model%ql(ip+1)) - model%glat(model%ql(ip)))
model%w(1:2) = model%w(1:2) * (1d0 - d)
model%w(3:4) = model%w(3:4) * d
end subroutine find_nearest

end program rads_add_era
