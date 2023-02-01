module tides_lpe
contains
!*lpetide -- Compute long-period equilibrium ocean tides
!+
subroutine lpetide (utc, lat, mode, lpe_tide, mf_tide)
use typesizes
use tides_aux
real(eightbytereal), intent(in) :: utc, lat
integer(fourbyteint), intent(in) :: mode
real(eightbytereal), intent(out) :: lpe_tide, mf_tide
!
! Computes the long-period equilibrium ocean tides.
! Fifteen tidal spectral lines from the Cartwright-Tayler-Edden
! tables or a total of 123 second and third order waves
! are summed over to compute the long-period tide.
!
! The terms with monthly, fortnightly, tri-monthly, weekly frequencies
! (Mm, Mf, Mtm and MSqm) and their sidelines are provided separately
! so they can be removed from models that have the combined
! equilibrium and non-equilibrium tides for these frequencies
! already in their solution (like FES2004)
!
! Furthermore, one can select to use only the 15 CTE
! waves or all 123 waves.
! mode = 0 : Use the 15 Cartwright-Taylor-Edden waves
!      = 1 : Use all 123 second and third order waves
!
! Technical references:
!   Cartwright & Tayler, Geophys. J. R.A.S., 23, 45, 1971.
!   Cartwright & Edden, Geophys. J. R.A.S., 33, 253, 1973.
!   Tamura Y., Bull. d'information des marées terrestres, Vol. 99, 1987.
!
! Based on LPEQMT by Richard Ray.
!
! Input arguments:
!  utc      : UTC time (seconds since 1985)
!  lat      : Latitude (degrees)
!  mode     : Determines how many waves are include (see above)
!
! Returned value:
!  lpe_tide : Computed long-period tide (meters)
!  mf_tide  : Only Mm, Mf, Mtm and MSqm tides and their sidelines
!             (meters)
!
!-
! $Log: lpetide.f90,v $
! Revision 1.9  2017/10/12 09:17:55  rads
! - Add generic mean_longtitudes routine
!
! Revision 1.8  2011/05/24 14:48:34  rads
! - Created module "tides"
!
! Revision 1.7  2009/05/26 00:18:01  rads
! - Include also Mf"' (2s-3n) to the mf_tide output. Minuscule effect on result.
!
! Revision 1.6  2009-05-20 11:51:46  rads
! - Add 8 more sidelines to mf_tide (for a total of 12)
!
! Revision 1.5  2008-06-03 17:32:32  rads
! - Changed from function to subroutine. Separately return Mm,Mf,Mtm,MSqm.
!
! Revision 1.4  2008-05-21 12:21:11  rads
! - Astronomical angles update to J2000 (same as gottide.f90)
!
! Revision 1.3  2008-03-16 09:55:06  rads
! - Ported from Fortran77 to Fortran90
!
! Revision 1.2  2006/09/25 00:22:38  rads
! - New code taken from FES2004 (Oct 2006 version)
! - Now uses 123 second and third order waves
!
! Revision 1.1  2005/10/27 15:23:43  rads
! - lpetide introduced as replacement for LPEQMT
!-----------------------------------------------------------------------
real(eightbytereal), parameter :: pi=3.14159265358979d0, rad=pi/180, &
	h2=0.609d0, k2=0.302d0, h3=0.291d0, k3=0.093d0, &
	h2k2=(1d0-h2+k2)*sqrt(5d0/4d0/pi)/1d5, h3k3=(1d0-h3+k3)*sqrt(7d0/4d0/pi)/1d5
integer(fourbyteint), parameter :: nw2=106,nw3=17
real(eightbytereal) :: shpn(5),sinlat,v
integer(fourbyteint) :: i,i0,i1

! Define the parameters for the second order waves
! s,   h,   p,  n',  p',    pot
integer(fourbyteint), parameter :: w2(6,nw2)=reshape((/ &
  4,  -2,   0,   0,   0,   -205, & ! MSqm
  4,  -2,   0,   1,   0,    -85, & ! MSqm'
  4,  -2,   0,   2,   0,     -8, & ! MSqm"
  3,   0,  -1,   2,   0,    -51, & ! Mtm"
  2,   0,   0,   3,   0,      7, & ! Mf"'
  3,   0,  -1,   0,   0,  -1275, & ! Mtm    ! CTE
  3,   0,  -1,   1,   0,   -528, & ! Mtm'   ! CTE
  2,   0,   0,   0,   0,  -6662, & ! Mf     ! CTE
  2,   0,   0,   1,   0,  -2762, & ! Mf'    ! CTE
  2,   0,   0,   2,   0,   -258, & ! Mf"    ! CTE
  1,   0,  -1,   0,   0,  -3518, & ! Mm     ! CTE
  1,   0,  -1,  -1,   0,    231, & ! Mm'    ! CTE
  1,   0,  -1,   1,   0,    228, & ! Mm"    ! CTE
  0,   0,   0,   1,   0,   2793, & ! CTE
  0,   1,   0,   0,  -1,   -492, & ! CTE
  0,   2,   0,   0,   0,  -3095, & ! CTE
  1,  -2,   1,   0,   0,   -673, & ! CTE
  2,  -2,   0,   0,   0,   -582, & ! CTE
  2,   0,  -2,   0,   0,   -288, & ! CTE
  3,  -2,   1,   0,   0,   -242, & ! CTE
  0,   0,   0,   2,   0,    -27, &
  0,   0,   2,   1,   0,      4, &
  0,   1,   0,  -1,  -1,     -4, &
  0,   1,   0,   0,   1,     26, &
  0,   1,   0,   1,  -1,      5, &
  0,   2,  -2,  -1,   0,      2, &
  0,   2,  -2,   0,   0,    -31, &
  0,   2,   0,   0,  -2,     -8, &
  0,   2,   0,   1,   0,     77, &
  0,   2,   0,   2,   0,     17, &
  0,   3,   0,   0,  -1,   -181, &
  0,   3,   0,   1,  -1,      3, &
  0,   4,   0,   0,  -2,     -7, &
  1,  -3,   1,  -1,   1,      2, &
  1,  -3,   1,   0,   1,    -29, &
  1,  -3,   1,   1,   1,      2, &
  1,  -2,  -1,  -2,   0,      3, &
  1,  -2,  -1,  -1,   0,      7, &
  1,  -2,   1,  -1,   0,     48, &
  1,  -2,   1,   1,   0,     43, &
  1,  -1,  -1,  -1,   1,      2, &
  1,  -1,  -1,   0,   1,    -21, &
  1,  -1,  -1,   1,   1,      0, &
  1,  -1,   0,   0,   0,     20, &
  1,  -1,   1,   0,  -1,      5, &
  1,   0,  -1,  -2,   0,     -3, &
  1,   0,   1,   0,   0,    189, &
  1,   0,   1,   1,   0,     77, &
  1,   0,   1,   2,   0,     21, &
  1,   1,  -1,   0,  -1,     18, &
  1,   2,  -1,   0,   0,     49, &
  1,   2,  -1,   1,   0,     24, &
  1,   2,  -1,   2,   0,      4, &
  1,   3,  -1,   0,  -1,      3, &
  2,  -4,   2,   0,   0,    -11, &
  2,  -3,   0,   0,   1,    -38, &
  2,  -3,   0,   1,   1,      2, &
  2,  -2,   0,  -1,   0,    -42, &
  2,  -2,   0,   1,   0,     37, &
  2,  -2,   2,   0,   0,      4, &
  2,  -1,  -2,   0,   1,     -4, &
  2,  -1,  -1,   0,   0,      3, &
  2,  -1,   0,   0,  -1,      7, &
  2,  -1,   0,   0,   1,    -20, &
  2,  -1,   0,   1,   1,     -4, &
  2,   0,  -2,  -1,   0,     15, &
  2,   0,  -2,   1,   0,     19, &
  2,   1,  -2,   0,  -1,      3, &
  2,   1,   0,   0,  -1,     23, &
  2,   1,   0,   1,  -1,      6, &
  2,   2,  -2,   0,   0,     20, &
  2,   2,  -2,   1,   0,      8, &
  2,   2,   0,   2,   0,      3, &
  3,  -5,   1,   0,   1,     -2, &
  3,  -4,   1,   0,   0,    -17, &
  3,  -3,  -1,   0,   1,     -7, &
  3,  -3,   1,   0,   1,    -12, &
  3,  -3,   1,   1,   1,     -4, &
  3,  -2,  -1,  -1,   0,    -10, &
  3,  -2,  -1,   0,   0,    -91, &
  3,  -2,  -1,   1,   0,      6, &
  3,  -2,   1,   1,   0,   -100, &
  3,  -2,   1,   2,   0,     -9, &
  3,  -1,  -1,   0,   1,    -13, &
  3,  -1,  -1,   1,   1,     -4, &
  3,  -1,   0,   0,   0,      6, &
  3,  -1,   0,   1,   0,      3, &
  3,  -1,   1,   0,  -1,      3, &
  3,   0,  -3,   0,   0,    -23, &
  3,   0,  -3,   1,  -1,      4, &
  3,   0,  -3,   1,   1,      4, &
  3,   0,   1,   2,   0,      5, &
  3,   0,   1,   3,   0,      2, &
  3,   1,  -1,   0,  -1,     11, &
  3,   1,  -1,   1,  -1,      4, &
  4,  -4,   0,   0,   0,     -8, &
  4,  -4,   2,   0,   0,     -6, &
  4,  -4,   2,   1,   0,     -2, &
  4,  -3,   0,   0,   1,    -14, &
  4,  -3,   0,   1,   1,     -6, &
  4,  -2,  -2,   0,   0,    -11, &
  4,  -1,  -2,   0,   1,     -3, &
  4,  -1,   0,   0,  -1,      3, &
  4,   0,  -2,   0,   0,   -169, &
  4,   0,  -2,   1,   0,    -70, &
  4,   0,  -2,   2,   0,     -6 /),(/6,nw2/))

! Define the parameters for the third order waves
! s,   h,   p,  n',  p',    pot
integer(fourbyteint), parameter :: w3(6,nw3)=reshape((/ &
  0,   0,   1,   0,   0,    -21, &
  0,   2,  -1,   0,   0,     -4, &
  1,  -2,   0,   0,   0,      4, &
  1,   0,   0,  -1,   0,     19, &
  1,   0,   0,   0,   0,   -375, &
  1,   0,   0,   1,   0,    -59, &
  1,   0,   0,   2,   0,      5, &
  2,  -2,   1,   0,   0,    -12, &
  2,   0,  -1,   0,   0,    -61, &
  2,   0,  -1,   1,   0,    -10, &
  3,  -2,   0,   0,   0,    -10, &
  3,   0,  -2,   0,   0,     -7, &
  3,   0,   0,   0,   0,    -30, &
  3,   0,   0,   1,   0,    -19, &
  3,   0,   0,   2,   0,     -4, &
  4,   0,  -1,   0,   0,     -8, &
  4,   0,  -1,   1,   0,     -5 /),(/6,nw3/))

! Compute 5 principal mean longitudes in radians at time t

call mean_longitudes (utc/86400d0, shpn(1), shpn(2), shpn(3), shpn(4), shpn(5))
shpn = shpn * rad

! Assemble long-period tide potential from the selected portion of waves
! Nodal term is included but not the constant term.
! Compute the potential of the second order waves

select case (mode)
case (0)	! CTE only
   i0=6 ; i1=20
case default	! All waves
   i0=1 ; i1=nw2
end select

! Compute the potential of the second order waves

sinlat = sin(lat*rad)

v = 0d0
do i=i0,13  ! First the waves that are already in FES2004
    v = v + w2(6,i) * cos(dot_product(w2(1:5,i),shpn))
enddo
mf_tide = h2k2 * v * (1.5d0*sinlat*sinlat - 0.5d0)

do i=14,i1  ! Then the rest
    v = v + w2(6,i) * cos(dot_product(w2(1:5,i),shpn))
enddo
lpe_tide = h2k2 * v * (1.5d0*sinlat*sinlat - 0.5d0)

if (mode == 0) return

! Compute the potential of the third order waves

v = 0d0
do i=1,nw3
    v = v + w3(6,i) * sin(dot_product(w3(1:5,i),shpn))
enddo

lpe_tide = lpe_tide + h3k3 * v * (2.5d0*sinlat*sinlat - 1.5d0) * sinlat
end subroutine lpetide

end module tides_lpe
