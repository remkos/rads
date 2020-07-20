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

!*rads_add_ecmwf -- Add ECMWF meteo models to RADS data
!+
! This program adjusts the contents of RADS altimeter data files
! with values computed from ECMWF operational meteorological models.
! The models provide sea level pressure, wet tropospheric correction, or
! wind speed.
!
! Input grids are found in the directory ${ALTIM}/data/ecmwf.
!
! Interpolation is performed in 6-hourly reduced Gaussian grids (N640);
! bi-linear in space, linear in time; contained in 6-hourly files.
!
! usage: rads_add_ecmwf [data-selectors] [options]
!
! References:
!
! Saastamoinen, J. (1972), Atmospheric corrections for the troposphere and
! stratosphere in radio ranging of satellites, in The Use of Artificial
! Satellites for Geodesy, Geophys. Monogr. Ser., vol. 15, edited by
! S. W. Hendriksen, A. Mancini, and B. H. Chovitz, pp. 247-251,
! American Geophysical Union, Washington, D.C.
!
! Petit, G., and B. Luzum (Eds.) (2010), IERS Conventions (2010),
! IERS Technical Note 36, Verlag des Bundesamts für Kartographie und Geodäsie.
!-----------------------------------------------------------------------
program rads_add_ecmwf

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
logical :: dry_on=.false., wet_on=.false., ib_on=.false., air_on=.false., wind_on=.false., new=.false., error

! Model data

integer(fourbyteint), parameter :: np_max=2140702, nj_max=1280
type :: model_
	integer(fourbytereal) :: nj, np, ip, pl(nj_max), ql(nj_max), idx(4)
	real(eightbytereal) :: glat(np_max), glon(np_max), w(4)
	real(eightbytereal) :: slp(np_max), wet(np_max), u10(np_max), v10(np_max)
end type
type(model_) :: m1, m2

! Initialise

call synopsis ('--head')
call rads_set_options ('dwiaun dry wet air ib all wind new')
call rads_init (S)

! Get ${ALTIM}/data/ecmwf/ directory

call parseenv ('${ALTIM}/data/ecmwf/', path)

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
	case ('u', 'wind')
		wind_on = .true.
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
if (rads_version ('Add ECMWF meteo models to RADS data', flag=flag)) return
call synopsis_devel (' [processing_options]')
write (*,1310)
1310  format (/ &
'Additional [processing_options] are:'/ &
'  -d, --dry                 Add ECMWF dry tropospheric correction' / &
'  -w, --wet                 Add ECMWF wet tropospheric correction' / &
'  -a, --air                 Add air tide' / &
'  -i, --ib                  Add static inverse barometer correction' / &
'  --all                     All of the above' / &
'  -u, --wind                Add ECMWF wind speed' / &
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
	u10(n), v10(n), f1, f2, g1, g2, slp, dslp, slp0

call log_pass (P)

! If "new" option is used, write only when fields are not yet available

ncid = P%fileinfo(1)%ncid
if (new .and. nff(nf90_inq_varid(ncid,'dry_tropo_ecmwf',i)) .and. &
	nff(nf90_inq_varid(ncid,'wet_tropo_ecmwf',i))) then
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

! Correct sea level pressure grids for altitude over land and lakes using standard atmosphere
! Convert pressure (mbar) to simple IB (m) if requested and over ocean
! Convert sea level pressure (mbar) to dry correction (m) using Saastamoinen models as referenced
! in IERS Conventions Chap 9

	if (surface_type(i) > 1.5d0) then ! Land and lakes
		g1 = -22768d-7 / (1d0 - 266d-5*cos(lat(i)*rad2) - 28d-8*h(i)) * pp_isa(h(i))
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

! Interpolate wet tropospheric correction in space and time

	if (wet_on) then
		wet(i) = f1 * dot_product (m1%w,m1%wet(m1%idx)) + f2 * dot_product (m2%w,m2%wet(m2%idx))
	endif

! Interpolate wind speed (both components) in space and time

	if (wind_on) then
		u10(i) = f1 * dot_product (m1%w,m1%u10(m1%idx)) + f2 * dot_product (m2%w,m2%u10(m2%idx))
		v10(i) = f1 * dot_product (m1%w,m1%v10(m1%idx)) + f2 * dot_product (m2%w,m2%v10(m2%idx))
	endif
enddo

! If no more fields are determined, abort.

if (.not.(dry_on .or. ib_on .or. wet_on)) then
	call log_records (0)
	stop
endif

! Store all data fields

call rads_put_history (S, P)

if (dry_on) call rads_def_var (S, P, 'dry_tropo_ecmwf')
if (wet_on) call rads_def_var (S, P, 'wet_tropo_ecmwf')
if (ib_on ) call rads_def_var (S, P, 'inv_bar_static')
if (air_on) call rads_def_var (S, P, 'dry_tropo_airtide')
if (wind_on) then
	call rads_def_var (S, P, 'wind_speed_ecmwf_u')
	call rads_def_var (S, P, 'wind_speed_ecmwf_v')
endif

if (dry_on) call rads_put_var (S, P, 'dry_tropo_ecmwf', dry)
if (wet_on) call rads_put_var (S, P, 'wet_tropo_ecmwf', wet)
if (ib_on ) call rads_put_var (S, P, 'inv_bar_static', ib)
if (air_on) call rads_put_var (S, P, 'dry_tropo_airtide', air)
if (wind_on) then
	call rads_put_var (S, P, 'wind_speed_ecmwf_u', u10)
	call rads_put_var (S, P, 'wind_speed_ecmwf_v', v10)
endif

call log_records (n)
end subroutine process_pass

!-----------------------------------------------------------------------
! get_gribs -- Load necessary ECMWF meteo grids
!-----------------------------------------------------------------------

function get_gribs (hex, model)
integer(fourbyteint), intent(in) :: hex
type(model_), intent(inout) :: model
logical :: get_gribs
!
! Input are mean sea level pressure files in the form
! ${ALTIM}/data/ecmwf/2016/msl_20160107_120000.grb
! and wet tropospheric correction files of the form
! ${ALTIM}/data/ecmwf/2016/wet_20160107_120000.grb
! and vector wind speed files of the form
! ${ALTIM}/data/ecmwf/2016/u10_20160107_120000.grb
! ${ALTIM}/data/ecmwf/2016/v10_20160107_120000.grb
!
! <hex> specifies the number of 6-hourly blocks since 1 Jan 1985.
! Data is stored in a buffer <model>
!-----------------------------------------------------------------------
logical :: get_lat

get_gribs = .true.
get_lat = .true.
if (dry_on .or. ib_on) then
	if (get_grib(hex,model,'msl',get_lat)) return
	get_lat = .false.
endif
if (wet_on) then
	if (get_grib(hex,model,'wet',get_lat)) return
	get_lat = .false.
endif
if (wind_on) then
	if (get_grib(hex,model,'u10',get_lat)) return
	get_lat = .false.
	if (get_grib(hex,model,'v10',get_lat)) return
endif
get_gribs = .false.
end function get_gribs

!-----------------------------------------------------------------------

function get_grib (hex, model, type, get_lat)
integer(fourbyteint), intent(in) :: hex
type(model_), intent(inout) :: model
character(len=3), intent(in) :: type
logical, intent(in) :: get_lat
logical :: get_grib
character(len=rads_cmdl) :: filenm
integer(fourbyteint) :: fileid, gribid, i, l, status, strf1985

600 format ('(',a,')')
1300 format (a,': ',a)

get_grib = .true.

l = strf1985(filenm, trim(path)//'%Y/'//type//'_%Y%m%d_%H.grb', hex*21600)
write (*,600,advance='no') filenm(l-18:l)

call grib_open_file(fileid,filenm,'r',status)
if (status /= grib_success) then
	write (*,1300) 'Error opening file',filenm(:l)
	return
endif
call grib_new_from_file(fileid,gribid,status)
if (status /= grib_success) then
	write (*,1300) 'Error loading grid',filenm(:l)
	return
endif

! Check sizes
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

! Read latitude band scructure
if (get_lat) then
	call grib_get(gribid,'pl',model%pl)
	model%ql(1) = 1
	do i = 2,model%nj
		model%ql(i) = model%ql(i-1) + model%pl(i-1)
	enddo
endif

select case (type)
case ('msl')
	call grib_get_data(gribid,model%glat,model%glon,model%slp)
case ('wet')
	call grib_get_data(gribid,model%glat,model%glon,model%wet)
case ('u10')
	call grib_get_data(gribid,model%glat,model%glon,model%u10)
case ('v10')
	call grib_get_data(gribid,model%glat,model%glon,model%v10)
end select

! Close file
call grib_release(gribid)
call grib_close_file(fileid)

! Successful completion
get_grib = .false.
end function get_grib

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

end program rads_add_ecmwf
