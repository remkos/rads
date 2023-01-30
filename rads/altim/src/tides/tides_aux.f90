module tides_aux
contains

!*mean_longitudes -- Compute 5 principle astronomical mean longitudes
!+
subroutine mean_longitudes (utc, s, h, p, np, pp)
use typesizes
real(eightbytereal), intent(in) :: utc
real(eightbytereal), intent(out) :: s, h, p, np, pp
!
! This routine computes the five principle mean longtitudes used to compute
! astronomical tides.
!
! Input is the UTC time in days since 1985. This is then converted to
! Julian centuries since 1.5 Jan 2000. Based on the mean longitudes and rates
! at that time, the mean longitudes at the current time are computed.
!
! Output are the mean longitudes of the moon, sun, lunar perigee, ascending
! lunar node, and solar perigee.
!
! Note: This routine uses time in UT and does not distinguish
! between the subtle differences of UTC, UT1, etc. However, this is
! more than accurate enough for our purposes.
! The formalae for the mean longitudes depend on dynamic time (DT).
! This routine assumes DT - UT = 63.48 seconds (2000).
!
! Reference:
! Jean Meeus, Astronomical Algorithms, 2nd ed., 1999.
!
! Arguments:
! utc     : UTC time in days since 1.0 January 1985
! s       : mean longitude of the moon (degrees)
! h       : mean longitude of the sun (degrees)
! p       : mean longitude of the lunar perigee (degrees)
! np      : negative of mean longitude of the ascending lunar node (degrees)
! pp      : mean longitude of the solar perigee (degrees)
!-
! $Log: tides_aux.f90,v $
! Revision 1.3  2019/03/08 09:57:29  rads
! - Introduced internal tide model HRET
!
! Revision 1.2  2017/12/04 15:03:32  rads
! - Update comments only
!
! Revision 1.1  2017/10/12 09:17:24  rads
! - Introduces tides_aux.f90 with generic tide routines, like mean_longitudes
!
! (c) Remko Scharroo - EUMETSAT
!-----------------------------------------------------------------------
real(eightbytereal) :: t

t  = (utc - 5478.4993d0) / 36525d0					! Julian centuries since 1.5 Jan 2000
s  = modulo(218.3164477d0 + 481267.88123421d0 * t, 360d0)	! Mean longitude of moon
h  = modulo(280.4662556d0 +  36000.76983081d0 * t, 360d0)	! Mean longitude of sun
p  = modulo( 83.3532465d0 +   4069.0137287d0  * t, 360d0)	! Mean longitude of lunar perigee
np = modulo(234.95548d0   +   1934.136261d0   * t, 360d0)	! Mean longitude of ascending lunar node
pp =        282.94d0      +      1.7192d0     * t		! Mean longitude of solar perigee

end subroutine mean_longitudes

end module tides_aux
