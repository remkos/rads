!****-------------------------------------------------------------------
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
!****-------------------------------------------------------------------

module rads_geo
use typesizes

real(eightbytereal), parameter :: ae_topex=6378136.3d0, f_topex=1d0/298.257d0
real(eightbytereal), parameter :: ae_wgs84=6378137.0d0, f_wgs84=1d0/298.257223563d0
real(eightbytereal), parameter :: ae_geosat=ae_wgs84  , f_geosat=f_topex

contains

!****f* rads_geo/checklon
! SUMMARY
! Convert and verify longitude
!
! SYNOPSIS
function checklon (lon_bounds, lon)
real(eightbytereal), intent(in) :: lon_bounds(2)
real(eightbytereal), intent(inout) :: lon
logical :: checklon
!
! PURPOSE
! This function converts the logitude <lon> (in degrees) to within limits
! <lon_bounds>, where lon_bounds(1) < lon_bounds(2). If lon is not within those
! boundaries, checklon is set to .true.
! The longitude will be wrapped (advanced or reduced by a multiple of 360
! degrees) if needed.
!
! If either of the limits is NaN, no wrapping and no check is performed.
!
! ARGUMENTS
! lon_bounds : Longitude limits (degrees)
! lon        : Longitude (degrees)
! checklon   : .true. if lon is outside limits, .false. otherwise.
!****-------------------------------------------------------------------
if (any(lon_bounds /= lon_bounds)) then
	checklon = .false.
	return
endif
lon = lon_bounds(1) + modulo (lon-lon_bounds(1),360d0)
checklon = (lon > lon_bounds(2))
end function checklon

!****f* rads_geo/sfdist
! SUMMARY
! Spherical distance between two points
!
! SYNOPSIS
function sfdist (lat0, lon0, lat1, lon1)
real(eightbytereal), intent(in) :: lat0, lon0, lat1, lon1
real(eightbytereal) :: sfdist
!
! PURPOSE
! Compute spherical distance between two points of given geocentric
! latitude and longitude. All angles are in radians.
!
! To avoid rounding errors at small distances, we use the haversine
! formula which is optimized for computing very small spherical
! distances, but only has inaccuracies determining distances to
! antipodal points.
!
! REFERENCES
! * R.W. Sinnott, "Virtues of the Haversine", Sky and Telescope,
!   vol. 68, no. 2, p. 159, 1984.
! * http://en.wikipedia.org/wiki/Great-circle_distance
!
! ARGUMENTS
! lat0, lon0 : Geocentic latitude and longitude of one point (rad)
! lat1, lon1 : Geocentic latitude and longitude of other point (rad)
! sfdist     : Spherical distance in radians.
!****-------------------------------------------------------------------
real(eightbytereal) :: s
s = sin((lat0-lat1)/2d0)**2 + cos(lat0)*cos(lat1)*sin((lon0-lon1)/2d0)**2
sfdist = 2d0*asin(sqrt(s))
end function sfdist

!****f* rads_geo/dh_ellipsoid
! SUMMARY
! Compute height difference between ellipsoids
!
! SYNOPSIS
function dhellips (conv, lat)
integer(fourbyteint), intent(in) :: conv
real(eightbytereal), intent(in) :: lat
real(eightbytereal) :: dhellips
!
! PURPOSE
! Compute height difference between WGS84, GEOSAT and TOPEX ellipsoids.
! The input is <lat> in degrees, the returned function value is <dhellips>
! in metres (the result is always non-negative). <conv> specifies the
! conversion (see below).
!
! The ellipsoids are defined by their equatorial radius and inverse
! flattening as shown in the table below. Also indicated are which
! products use which ellipsoids.
!
! Ellipsoid  Eq. radius (m)  Inv. flattening (-)  Products
! WGS84      6378137.0       298.257223563        ESA, Sentinel-3, Jason-CS
! GEOSAT     6378137.0       298.257              GEOSAT T2 GDRs,
!                                                 DEOS ERS orbits
! TOPEX      6378136.3       298.257              T/P and Jason GDRs,
!                                                 GEOSAT J3 GDRs, RADS
!
! The conversion is linearised in sin(lat)**2, changes in latitude
! are ignored.
!
! ARGUMENTS:
! conv     : Conversion indicator
!                1 : Height of WGS84 over TOPEX ellipsoid
!                2 : Height of WGS84 over GEOSAT ellipsoid
!                3 : Height of GEOSAT over TOPEX ellipsoid
!            Other : Zero height difference is returned
! lat      : Latitude (degrees)
! dhellips : Height difference between ellipsoids (m)
!****-------------------------------------------------------------------
real(eightbytereal), parameter :: pi = 3.1415926535897932d0, rad = pi/180d0

select case (conv)
case (1) ! WGS84 - TOPEX
	dhellips = (ae_wgs84 - ae_topex) + &
			(ae_topex*f_topex - ae_wgs84*f_wgs84) * sin(lat*rad)**2
case (2) ! WGS84 - GEOSAT
	dhellips = ae_geosat * (f_geosat - f_wgs84) * sin(lat*rad)**2
case (3) ! GEOSAT - TOPEX
	dhellips = (ae_geosat - ae_topex) * (1d0 - f_geosat * sin(lat*rad)**2)
case default ! Other
	dhellips = 0d0
end select
end function dhellips

end module rads_geo
