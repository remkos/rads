module tides_pole
use typesizes

! Polar motion table

integer(fourbyteint), parameter :: mjd0 = 40587, mjd2 = 62502	! 1970-2030
integer(fourbyteint) :: mjd1
logical :: poletideinit_done = .false.
type poletab_
	real(eightbytereal) :: xpole, ypole, utc_ut1
end type
type(poletab_), allocatable :: poletab(:)

! Polar motion and pole tide constants:
! Amp = amplitude in meters when Xpole and Ypole are in arcsec.
! Xpole0 = Average X position of Pole in arcsec
! Ypole0 = Average Y position of Pole in arcsec
! (X = along Greenwich meridian (towards equator), Y = towards 90 West)

real(eightbytereal), parameter :: amp=-69.435d-3,rad=atan(1d0)/45d0

contains

!*poletide -- Compute pole tide elevation
!+
function poletide (mjd, lat, lon)
real(eightbytereal), intent(in) :: mjd, lat, lon
real(eightbytereal) :: poletide
!
! This function computes the pole tide elevation (poletide) for a given
! location (lat and lon) and epoch (mjd). The position of the pole
! is interpolated from the daily "IERS 08 C04" values contained in
! $ALTIM/data/tables/eopc04.62-now  and
! $ALTIM/data/tables/eopc04_extended.dat
!
! These solutions are consistent with ITRF 2008 and originate from
! http://hpiers.obspm.fr/eoppc/eop/eopc04/eopc04.62-now  and
! http://hpiers.obspm.fr/eoppc/series/prediction/eopc04_extended.dat
!
! Although the tables provide values from 1962 till present, only
! values from 1970 onward are stored (currently limited to 2030).
! The daily values are interpolated linearly.
!
! The pole tide elevation algorithm is based on the references below.
! A linear drift was removed to take out the effect of GIA, which
! should be corrected for separately (Wahr et al. 2015, Desai, 2015).
!
! The pole tide provided is that of the solid earth plus ocean. This
! means that the value provided is (1 + k2) times the equilibrium pole
! tide (He) with k2 = 0.302. The respective contributions from the solid
! earth and ocean are h2 He (with h2 = 0.609) and (1 - h2 + k2) He,
! respectively. In other words, of the total 46.8% is due to the solid
! earth and 53.2% to the ocean.
!
! References:
! Munk, W. H., and G. J. F. MacDonald,
! The Rotation of the Earth: A Geophysical Discussion,
! Cambridge University Press, New York, 1960.
!
! Wahr, J. W., Deformation of the Earth induced by polar motion,
! J. Geophys. Res., 90, 9363-9368, 1985.
!
! Desai, S. D., Observing the pole tide with satellite altimetry,
! J. Geophys. Res., 107(C11), 3186, doi: 10.1029/2001JC001224, 2002.
!
! Bizouard, C., and D. Gambis, The combined solution C04 for Earth
! Orientation Parameters consistent with International Terrestrial
! Reference Frame 2008, Tech. rep., Observatoire de Paris, 2011.
!
! Wahr, J. M., R. S. Nerem, and S. V. Bettadpur, The pole tide and its
! effect on GRACE time-variable gravity measurements: implications for
! the estimates of surface mass variations, J. Geophys. Res., 120,
! doi:10.1002/2015JB011986, 2015.
!
! Desai, S. D., personal communication, 2015.
!
! Input arguments:
!  mjd      : Modified Julian Date (floating point double)
!  lat, lon : Geodetic latitude and longitude of the observer in degrees.
!
! Returned value:
!  poletide : Pole tide elevation (solid + ocean) in meters.
!-
! $Log: poletide.f90,v $
! Revision 1.18  2015/08/10 08:23:46  rads
! - Period increased to 1970-2030
!
! Revision 1.17  2015/08/10 08:15:07  rads
! - Removed the code to load the mean pole table
! - Introduced removal of linear trend in mean pole
!
! Revision 1.16  2015/06/25 14:07:16  rads
! - Update of comments
! - Add code to allow interpolation of mean pole table (not used yet)
!
! Revision 1.15  2014/10/17 07:35:56  rads
! - Fixed bug in v1.14 which skipped every second line
!
! Revision 1.14  2014/09/26 20:37:26  rads
! - Avoid failure to read file when number of header records changes.
! - Made to work until 2030 (hopefully)
!
! Revision 1.13  2012/02/14 01:51:55  rads
! - Fixed error causing segmentation fault
!
! Revision 1.12  2012/02/01 20:06:09  rads
! - Introduced new poletideinit routine
!
! Revision 1.11  2011/05/24 14:48:34  rads
! - Created module "tides"
!
! Revision 1.10  2009/05/06 21:02:57  rads
! - Use new prediction file (eopc04_extended.dat)
! - Let values in eopc04.62-now prevail over predictions
!
! Revision 1.9  2008/03/16 13:53:09  rads
! - Ported from Fortran77 to Fortran90
!
! Revision 1.8  2007/10/19 21:39:26  rads
! - Changed to use new EOPC04_05 tables (consistent with ITRF2005)
!
! Revision 1.7  2005/11/17 16:56:23  rads
! - Read directly from tables with daily values available from IERS
!
! 23-Feb-2004 - Extended pole tide table and included check
! 10-Apr-2000 - Added reference.
! 11-Feb-1997 - Created by Remko Scharroo
!-----------------------------------------------------------------------
integer(fourbyteint) :: i
real(eightbytereal) :: xpole,ypole,d

if (mjd1 < mjd0) stop "poletide: use poletideinit first"

! Interpolate the pole position table and compute the pole tide

i = int(mjd)
d = mjd - i
if (i < mjd0 .or. i >= mjd1) stop "poletide: epoch requested outside table"
xpole = poletab(i)%xpole * (1-d) + poletab(i+1)%xpole * d
ypole = poletab(i)%ypole * (1-d) + poletab(i+1)%ypole * d

poletide = amp * sin(2d0*lat*rad)*(xpole*cos(lon*rad)-ypole*sin(lon*rad))
end function poletide

subroutine poletideinit ()
integer(fourbyteint) :: i, unit, freeunit
character(128) :: filenm

if (poletideinit_done) stop "poletideinit: already allocated"

! Allocate memory for pole table and initialise

allocate (poletab(mjd0:mjd2))
poletab = poletab_(0d0,0d0,0d0)
mjd1 = mjd0 - 1 ! Just before allowed start of table

! Load two pole tables

unit = freeunit()
filenm = '/user/altim'
call checkenv('ALTIM',filenm,i)
call poletide_load (filenm(:i)//'/data/tables/eopc04.62-now')
call poletide_load (filenm(:i)//'/data/tables/eopc04_extended.dat')

poletideinit_done = .true.

contains

! Load the pole position table, after mjd1 and before mjd2

subroutine poletide_load (filenm)
character(*), intent(in) :: filenm
integer(fourbyteint) :: i,ios
real(eightbytereal) :: x,y,u,t
character(len=155) :: line
integer(fourbyteint), parameter :: mjd2000 = 51544

write (*,600) trim(filenm)
open (unit,file=filenm,status='old')
do
	read (unit,'(a)',iostat=ios) line
	if (ios /= 0) exit	! Stop reading at end of file
	if (line(:2) /= '19' .and. line(:2) /= '20') cycle ! Skip lines NOT starting with 19 or 20
	read (line,610,iostat=ios) i,x,y,u
	if (ios /= 0) cycle ! Skip erroneous lines
	if (i <= mjd1) cycle ! Skip lines already loaded (or before mjd0)
	if (i > mjd2) stop "poletide: too many table values"
	t = (i - mjd2000) / 365.25d0	! Years since 2000.0
	! Remove GIA contribution to polar motion (Desai, 2015)
	poletab(i)%xpole = x - (0.05097d0 + 0.00062d0 * t)
	poletab(i)%ypole = y - (0.33449d0 + 0.00348d0 * t)
	! Also store UTC-UT1 (not used)
	poletab(i)%utc_ut1 = u
	mjd1 = i
enddo
close (unit)
write (*,601)

600 format ('(poletide: loading ',a,$)
601 format (')')
610 format (12x,i7,3f11.6)
end subroutine poletide_load

end subroutine poletideinit

subroutine poletidefree ()
if (allocated(poletab)) deallocate(poletab)
poletideinit_done = .false.
end subroutine poletidefree

end module tides_pole
