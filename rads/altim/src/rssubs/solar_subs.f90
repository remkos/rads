module solar_subs

contains

subroutine sun_orbit_geometry (mjd, raan, incl, alpha_prime, beta_prime)
use typesizes
real(eightbytereal), intent(in) :: mjd, raan, incl
real(eightbytereal), intent(out) :: alpha_prime, beta_prime
!
! This routine computes the alpha' and beta' angles of the sun to any orbit
! with a given ! right ascension of the ascending node (raan) and
! inclination (incl).
! Alpha' is the angle between the ascending node and the 6:00 LST position
! of the satellite, i.e., where the earth-sun vector is perpendicular to
! the satellite radius and the satellite is flying into the sun.
! Beta' is the angle between the sun and the orbital plane.
!
! Input:
!  mjd   : UTC time in Modified Julian Dates
!  raan  : Right ascension of the ascending node (radian)
!  incl  : Inclination (radian)
!
! Output:
!  alpha_prime : alpha' angle (radian)
!  beta_prime  : beta' angle (radian)
!-
! 2010-04-14 - Created by Remko Scharroo (Altimetrics LLC)
!-----------------------------------------------------------------------
real(eightbytereal) :: n(3), s(3), a(3), b(3), r, ra, decl
real(eightbytereal), parameter :: pi = 4d0 * atan(1d0), rad = pi / 45d0
!
! Determine the direction of the ascending node
a = (/ cos(raan), sin(raan), 0d0 /)
!
! Determine the direction of the normal to the orbital plane
r = sin(incl)
n = (/ r * a(2), -r * a(1), cos(incl) /)
!
! Determine the direction of the sun
call solar_coord (mjd, r, ra, decl)
r = cos(decl)
s = (/ r * cos(ra), r * sin(ra), sin(decl) /)
!
! Compute the angle between the orbit and the sun
beta_prime = asin(dot_product(n, s))
!
! Compute the vector perpendicular to sun and orbit normal as (sun x normal) and normalise it.
! This is the orbit 6:00 LST location.
b = (/ s(2) * n(3) - s(3) * n(2), s(3) * n(1) - s(1) * n(3), s(1) * n(2) - s(2) * n(1) /)
b = b / sqrt(dot_product(b,b))
!
! Now compute alpha_prime
alpha_prime = sign(acos(dot_product(a,b)),b(3))
end subroutine sun_orbit_geometry

function gmst (mjd)
use typesizes
real(eightbytereal), intent(in) :: mjd
real(eightbytereal) :: gmst
!
! This function computes Greenwich mean sidereal time (expressed in radian)
! based on the terrestial time in MJD
!
! Source: http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-sider-time
!
! Input:
!  mjd   : UTC time in Modified Julian Dates
!
! Return value:
!  gmst  : Greenwich mean sidereal time in radian
!-
! 2010-04-14 -- Created by Remko Scharroo (Altimetrics LLC)
!-----------------------------------------------------------------------
real(eightbytereal) :: D
real(eightbytereal), parameter :: hr2rad = atan(1d0) / 3d0

! Days since J2000 in Terrestial Time
D = mjd - 51544.499257d0

! Compute GMST (changing angle in hours to radian)
gmst = modulo(18.697374558d0 + 24.06570982441908d0 * D, 24d0) * hr2rad
end function gmst

subroutine solar_coord (mjd, r, ra, decl)
use typesizes
real(eightbytereal), intent(in) :: mjd
real(eightbytereal), intent(out) :: r, ra, decl
!
! This routine computes the location of the sun in an earth centered rotating
! reference frame as a function of time.
!
! Source: http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-solar
!
! Input:
!  mjd   : UTC time in Modified Julian Dates
!
! Output:
!  r     : Distance to the sun in AU
!  ra    : Right-ascension of the sun in radians
!  decl  : Declination of the sun in radians
!-
! 2006-11-16 - Created by Remko Scharroo (Altimetrics LLC)
!-----------------------------------------------------------------------
real(eightbytereal), parameter :: rad = atan(1d0) / 45d0
real(eightbytereal) :: D, g, L, e, sinL

! Days since J2000 in Terrestial Time
D = mjd - 51544.499257d0
! Siderial time
g = (357.529d0 + 0.98560028d0 * D) * rad
! Ecliptic longitude corrected for aberration
L = (280.459d0 + 0.98564736d0 * D + 1.915d0 * sin (g) + 0.020d0 * sin (2*g)) * rad
! Sun-earth distance
r = 1.00014d0 - 0.01671d0 * cos (g) - 0.00014d0 * cos (2*g)
! Obliquity of the ecliptic
e = (23.439d0 - 36d-8 * d) * rad

sinL = sin(L)
! Right Ascension
ra = atan2 (cos(e) * sinL, cos(L))
! Declination
decl = asin (sin(e) * sinL)
end subroutine solar_coord

function eclipse (ra, decl, lon, lat, h)
use typesizes
real(eightbytereal), intent(in) :: ra, decl, lon, lat, h
logical :: eclipse
!
! This routine determines whether an object at location (lon, lat, h) is
! in the earth's shadow given the sun's position (ra, decl).
!
! Simplifications: spherical earth, cilindrical shadow
! So flattening of the earth, conical shape of the shadow and penumbra are
! ignored.
!
! Input:
!  ra, decl : Right-ascension and declination of the Sun (in radians)
!  lon, lat : Location's longitude and latitude (in radians)
!  h        : Locations altitude (in meters)
!
! Output:
!  eclipse  : TRUE if eclipsed, FALSE if in sunlight
!-
! 2006-11-16 - Created by Remko Scharroo (Altimetrics LLC)
!-----------------------------------------------------------------------
real(eightbytereal), parameter :: re = 6371.0d3
real(eightbytereal) :: alpha, g

! Half the angular width of the earth shadow at altitude h as seen from the
! earth centre
alpha = asin(1d0/(1d0+h/re))

! Compute dot product of coordinates of location and anti-sun.
g = -cos(decl) * cos(lon-ra) * cos(lat) - sin(decl) * sin(lat)

! If dot product > cos(alpha) then where in shadow.
eclipse = (g > cos(alpha))
end function eclipse

end module solar_subs
