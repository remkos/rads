!****-------------------------------------------------------------------
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
!****-------------------------------------------------------------------

module rads_geo
use typesizes

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
! lon_bounds : longitude limits (degrees)
! lon        : longitude (degrees)
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

end module rads_geo
