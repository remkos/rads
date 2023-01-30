module tides_fes
use typesizes
integer(twobyteint), parameter, private :: undef=-32768
integer(fourbyteint), parameter, private :: nw=33
real(eightbytereal), parameter, private :: pi=4d0*atan(1d0),rad=pi/180d0,xmin=0d0,xmax=360d0,ymin=-90d0,ymax=90d0

! These values are internal to the program and to the nodal correction
! routines. Do not change under any circumstances.
!  9 Primary tides
!  2 Additional small tides
!  4 Long-period tides
! 18 Extra diurnal and semi diurnal infered by admittance

character(8), parameter, private :: wave(nw) = (/ &
	'Q1      ','O1      ','K1      ','2N2     ','N2      ','M2      ','K2      ', &
	'S2      ','P1      ','M4      ','S1      ','Mf      ','Mm      ','Mtm     ', &
	'MSqm    ','Nu2     ','Mu2     ','L2      ','T2      ','Eps2    ','Lambda2 ', &
	'Eta2    ','2Q1     ','Sigma1  ','Ro1     ','M11     ','M12     ','Ki1     ', &
	'Pi1     ','Phi1    ','Teta1   ','J1      ','OO1     ' /)
real(eightbytereal), parameter, private :: freq(nw) = (/ 13.39866087990d0, 13.94303558000d0, 15.04106864000d0, 27.89535481990d0, &
	28.43972952010d0, 28.98410422000d0, 30.08213728000d0, 30.00000000000d0, 14.95893136000d0, &
	57.9682084d0, 15.0000000d0, 1.09803310d0, 0.54437470d0, 1.64240780d0, 2.11392880d0, &
	28.51258314000d0, 27.96820844000d0, 29.52847892000d0, 29.95893332010d0, 27.4238337d0, 29.4556253d0, &
	30.6265120d0, 12.8542862d0, 12.9271398d0, 13.4715145d0, 14.4966939d0, 14.4874103d0, &
	14.5695476d0, 14.9178647d0, 15.1232059d0, 15.5125897d0, 15.5854433d0, 16.1391017d0 /) * rad

type :: festideinfo
	integer(fourbyteint) ::	nx,ny,ng
	logical :: haveload
	integer(twobyteint), allocatable :: p_ocean(:,:,:,:),p_load(:,:,:,:)
	real(eightbytereal) :: dx,dy,z0(nw,2),dz(nw,2),f(nw),v0_u(nw),t_nodal,nan,pmin
endtype

contains

!*festide -- Compute tides according to FES model
!+
subroutine festide (info, utc, lat, lon, otide_sp, otide_lp, ltide_sp, ltide_lp)
type(festideinfo), intent(inout) :: info
real(eightbytereal), intent(in) :: utc, lat, lon
real(eightbytereal), intent(out) :: otide_sp, otide_lp, ltide_sp, ltide_lp

! This routine makes the tidal predictions of ocean and load tide
! (optional) based on one of the FES models. This routine is heavily
! based on the routines by J.M. Molines and F. Lefevre for the FES99,
! FES2002 and FES2004 models.
!
! The input grids can be found in $ALTIM/data/<tide>. The grids are
! in NetCDF format and were converted using the program 'fes2nc'.
!
! To initialize the computation, the subroutine festideinit should be
! called first. It allocates the appropriate amount of memory and
! loads the grids into memory. To release the memory for further
! use, call festidefree.
!
! Longitude and latitude are to be specified in degrees; time in UTC
! seconds since 1 Jan 1985. All predicted tides are output in meters.
! If the tide is requested in a point where it is not defined, NaN
! (Not-a-Number) is returned.
!
! Note that the long-period tides (tide_lp) include only the
! waves (maximum 4) that are included in the FES model. Use the
! routine lpetide to obtain all long-period waves.
!
! Input arguments:
!  info     : Struct initialised by festideinit
!  utc      : UTC time in seconds since 1 Jan 1985
!  lat      : Latitude (degrees)
!  lon      : Longitude (degrees)
!
! Output arguments:
!  otide_sp : Predicted short-period ocean tide (m)
!  otide_lp : Predicted long-period ocean tide (m)
!  ltide_sp : Predicted short-period load tide (m)
!  ltide_lp : Predicted long-period load tide (m)
!-
! $Log: festide.f90,v $
! Revision 1.34  2019/03/08 08:36:23  rads
! - Making parameters private
!
! Revision 1.33  2014/09/12 13:51:06  rads
! - Added option to use FES2012
!
! Revision 1.32  2014/05/03 16:21:46  rads
! - Allow changing of interpolation limit using info%pmin
!
! Revision 1.31  2014/02/08 16:39:01  rads
! - Replaced cmplx by dcmplx to avoid loss of precision
!
! Revision 1.30  2011/05/24 14:47:48  rads
! - Created module "tides"
! - Assign "info" argument on initialisation
! - Renamed *free to *tidefree, *init to *tideinit
!
! Revision 1.29  2009/05/20 21:10:03  rads
! - Added additional argument: ltide_lp
!
! Revision 1.28  2009/05/20 15:54:56  rads
! - FES2004 now uses 15 load grids
!
! Revision 1.27  2009/05/12 20:55:39  rads
! - All tidal components now stored in one grid
!
! Revision 1.26  2009/01/16 16:35:34  rads
! - Arrays of characters need to be initialized with same length (gfortran 4.4)
!
! Revision 1.25  2008/06/03 21:29:47  rads
! - No longer computes equilibrium tides.
!   Only returns long-period waves in model.
!
! Revision 1.24  2008/03/17 17:17:43  rads
! - Ported from Fortran 77 to 90. Speed performance improvements
!
! Revision 1.23  2007/04/02 17:30:50  rads
! - Update to 2007 version of FES2004
!
! Revision 1.22  2006/09/29 15:08:50  rads
! - Send progress reports to standard error, not standard output
!
! Revision 1.21  2006/09/25 01:51:45  rads
! - Use netCDF grids as input
! - Fix amplitude of M4 tide (was about 9% too small)
!
! Revision 1.18  2006/08/07 17:24:59  rads
! - Removing obsolete Fortran code so it can be used with gfortran
!
! Revision 1.15  2005/10/27 20:14:43  rads
! - Added S1 and M4 short period tides and Mm, Mf, Mtm and MSqm long period tides
!
! Revision 1.13  2005/05/31 20:52:45  rads
! - Do all internal computation in double precision
!
! Revision 1.12  2005/03/10 10:57:44  eelco
! Removed calls to perror and added makefiles for compatibility with Intel Fortran Compiler 8.1
!
! Revision 1.11  2004/11/23 01:54:11  remko
! - Implementation supports FES2004
! - Automatically allocates memory, using FESTIDEINIT and FESTIDEFREE
! - Change of arguments in FESTIDE
!
! 10-Jun-2003 - Prepared for FES2002 (P1 component from grids);
!               switched M11/M12 admittance amplitudes;
!               fixed M12 frequency; allow grids with various scales.
!  9-Oct-2001 - Adapted by Remko Scharroo for DEOS
! 09-May-2001 - FES99 version by Fabien Lefevre (CLS)
! 16-Jun-1995 - First version by J. M. Molines (Grenoble)
!-----------------------------------------------------------------------
real(eightbytereal), parameter :: delta_max=24d0
real(eightbytereal) :: slon,delta,t,fwv(nw)
complex(eightbytereal) :: val(nw),w(nw)
integer(fourbyteint) :: istat1,istat2

! When called with invalid time, reset time reference and return

if (utc > 1d20) then
	info%t_nodal=1d30
	return
endif

! If latitude out of range, bail out

if (lat < ymin .or. lat > ymax) then
	otide_sp = info%nan
	otide_lp = info%nan
	ltide_sp = info%nan
	ltide_lp = info%nan
	return
endif

! Limit longitude to interval of grids

slon=lon
if (slon > xmax) then
	slon = slon-360d0
else if (slon < xmin) then
	slon = slon+360d0
endif

! Conversion sec85 to days since 1900

t = utc/86400d0+31046d0

! Determine time since last call to fes_astron (in hours).
! If time lapse is larger than delta_max, call astronimics again.

delta = (t-info%t_nodal)*24d0
if (abs(delta) >= delta_max) then
	call fes_astron(t/36525d0)
	info%t_nodal = t
	delta = 0d0
endif

! Compute all ocean tide components

call fes_interp(info%p_ocean,info%dz(:,1),info%z0(:,1),slon,lat,val,istat1)
if (istat1 == 0) then
	otide_sp = info%nan
	otide_lp = info%nan
else
	w = expj(-(freq*delta+info%v0_u))
	fwv = info%f*dble(w*val)
	otide_sp = sum(fwv(1:11)) + sum(fwv(16:nw))
	otide_lp = sum(fwv(12:15))
endif

! Compute all load tide components (if requested)

if (.not.info%haveload) then
	ltide_sp = 0d0
	ltide_lp = 0d0
	return
endif

call fes_interp(info%p_load,info%dz(:,2),info%z0(:,2),slon,lat,val,istat2)
if (istat2 == 0) then
	ltide_sp = info%nan
	ltide_lp = info%nan
else
	if (istat1 == 0) w = expj(-(freq*delta+info%v0_u))
	fwv = info%f*dble(w*val)
	ltide_sp = sum(fwv(1:11)) + sum(fwv(16:nw))
	ltide_lp = sum(fwv(12:15))
endif

contains
!-----------------------------------------------------------------------
! Perform bi-linear interpolation at point x,y from the gridded files.
! istat returns the number of points used for the interpolation

subroutine fes_interp(work,fact,offs,x,y,val,istat)
integer(twobyteint), intent(in) :: work(:,:,:,:)
real(eightbytereal), intent(in) :: x,y,fact(:),offs(:)
complex(eightbytereal), intent(out) :: val(:)
integer(fourbyteint), intent(out) :: istat
integer(fourbyteint) :: kw,i0,j0,i,j
real(eightbytereal) :: ptot,xij,yij,pds

istat = 0

! Compute indices for desired position

xij = (x-xmin)/info%dx+1
yij = (y-ymin)/info%dy+1
i0 = min(info%nx-1,int(xij))
j0 = min(info%ny-1,int(yij))

! Sum up the weighted values for the corners

ptot = 0d0
val = (0d0,0d0)

do i = i0,i0+1
	jloop: do j = j0,j0+1
! Check if all major constituents are available:
		do kw=1,info%ng
			if (work(i,j,1,kw) == undef) cycle jloop
		enddo

		istat = istat+1
		pds = (1d0-abs(i-xij))*(1d0-abs(j-yij))
		ptot = ptot+pds

		forall (kw=1:info%ng) val(kw)=val(kw)+(work(i,j,1,kw)*fact(kw)+offs(kw))*expj(work(i,j,2,kw)*1d-2*rad)*pds
	enddo jloop
enddo

if (istat == 0 .or. ptot <= info%pmin) then
	istat=0
	return
endif
val(1:info%ng)=val(1:info%ng)/ptot

! Infer additional constituents by admittance

! SEMI-DIURNAL (from Grenoble to take advantage of 2N2)

val(16) = -0.0061047d0*val(7)+0.1568788d0*val(5)+0.0067557d0*val(6)	! Nu2 from N2,M2,K2
val(17) =  0.0694400d0*val(7)+0.3515356d0*val(5)-0.0462783d0*val(6)	! Mu2 from N2,M2,K2
val(18) =  0.0771378d0*val(7)-0.0516535d0*val(5)+0.0278699d0*val(6)	! L2 from N2,M2,K2
val(19) =  0.1804802d0*val(7)-0.0201012d0*val(5)+0.0083315d0*val(6)	! T2 from N2,M2,K2
val(20) =  0.53285d0  *val(4)-0.03304d0  *val(5)				! Eps2 from 2N2,N2
val(21) =  0.0165036d0*val(7)-0.0133078d0*val(5)+0.0077534d0*val(6)	! Lambda2 from N2,M2,K2
val(22) = -0.0034925d0*val(6)+0.0831707d0*val(7)				! Eta2 from M2,K2

! DIURNALS (from Richard Ray perth2 program)

val(23) =  0.263d0 *val(1)-0.0252d0*val(2)				! 2Q1 from Q1,O1
val(24) =  0.297d0 *val(1)-0.0264d0*val(2)				! Sigma1 from Q1,O1
val(25) =  0.164d0 *val(1)+0.0048d0*val(2)				! Ro1 from Q1,O1
val(26) =  0.0389d0*val(2)+0.0282d0*val(3)				! M11 from O1,K1
val(27) =  0.0140d0*val(2)+0.0101d0*val(3)				! M12 from O1,K1
val(28) =  0.0064d0*val(2)+0.0060d0*val(3)				! K1 from O1,K1
val(29) =  0.0030d0*val(2)+0.0171d0*val(3)				! Pi1 from O1,K1
val(30) = -0.0015d0*val(2)+0.0152d0*val(3)				! Phi1 from O1,K1
val(31) = -0.0065d0*val(2)+0.0155d0*val(3)				! Teta1 from O1,K1
val(32) = -0.0389d0*val(2)+0.0836d0*val(3)				! J1 from O1,K1
val(33) = -0.0431d0*val(2)+0.0613d0*val(3)				! OO1 from O1,K1

! P1 from Grenoble admittance code when grid is not provided

if (info%ng < 9) &
val( 9) = -0.2387302d0*val(1)+0.1038608d0*val(2)+0.2892755d0*val(3)	! P1 from Q1,O1,K1

end subroutine fes_interp

!-----------------------------------------------------------------------
! Initialize some astronomic data and compute nodal corrections.

subroutine fes_astron(tj)
real(eightbytereal), intent(in) :: tj
real(eightbytereal) :: n,p,s,p1,pp,nu,xi,tt,nuprim,nusec,hp,r,iang,x1ra,hpi,tgn2,at1,at2,u,tanhn,coshn,sinhn,sin1n,sin2n

tt = mod(180d0         +  360d0*36525d0*tj, 360d0)*rad	! Mean solar angle relative to Greenwich
n  = mod(259.1560563d0 - 1934.1423972d0*tj, 360d0)*rad	! Longitude of ascending lunar node
hp = mod(280.1895015d0 +  36000.76892d0*tj, 360d0)*rad	! Mean solar longitude
s  = mod(277.0256206d0 +   481267.892d0*tj, 360d0)*rad	! Mean lunar longitude
p1 = mod(281.2208568d0 +     1.719175d0*tj, 360d0)*rad	! Longitude of solar perigee
p  = mod(334.3837214d0 + 4069.0322056d0*tj, 360d0)*rad	! Longitude of lunar perigee

u = 9.13694997d-1-3.5692561d-2*cos(n)
iang = acos(u)
tgn2 = tan(n/2)
at1 = atan(1.01883d0*tgn2)
at2 = atan(6.4412d-1*tgn2)
xi = -at1-at2+n
if (n > pi) xi = xi-2*pi
nu = at1-at2

sinhn = sin(iang/2d0)
sin1n = sin(iang)
sin2n = sin(2d0*iang)
coshn = cos(iang/2d0)
tanhn = sinhn/coshn

! For constituents L2, K1, K2

pp=p-xi
x1ra = sqrt(1-12*tanhn**2*cos(2*pp)+36*tanhn**4)
r = atan(sin(2*pp)/(1/(6*tanhn**2)-cos(2*pp)))
nuprim = atan(sin2n*sin(nu)/(sin2n*cos(nu)+3.347d-1))
nusec = 0.5d0*atan(sin1n**2*sin(2*nu)/(sin1n**2*cos(2*nu)+7.27d-2))

! Compute nodal corrections from Schureman (1958)

info%f( 1) = sin1n*coshn**2/0.38d0
info%f( 2) = info%f( 1)
info%f( 3) = sqrt(0.8965d0*sin2n**2+0.6001d0*sin2n*cos(nu)+0.1006d0)
info%f( 4) = coshn**4/0.9154d0
info%f( 5) = info%f( 4)
info%f( 6) = info%f( 4)
info%f( 7) = sqrt(19.0444d0*sin1n**4+2.7702d0*sin1n**2*cos(2*nu)+.0981d0)
info%f( 8) = 1d0
info%f( 9) = info%f( 8)
info%f(10) = info%f( 4)**2
info%f(11) = 1d0
info%f(12) = sin1n**2/0.1578d0
info%f(13) = (2d0/3d0-sin1n**2)/0.5021d0
info%f(14) = info%f(12)
info%f(15) = info%f(12)
info%f(16) = info%f( 4)
info%f(17) = info%f( 4)
info%f(18) = info%f( 4)*x1ra
info%f(19) = info%f( 8)
info%f(20) = info%f( 8)
info%f(21) = info%f( 4)
info%f(22) = sin1n**2/0.1565d0
info%f(23) = info%f( 1)
info%f(24) = info%f( 1)
info%f(25) = info%f( 1)
info%f(26) = sin2n/0.7214d0
info%f(27) = info%f( 1)
info%f(28) = info%f(26)
info%f(29) = info%f( 8)
info%f(30) = info%f( 8)
info%f(31) = info%f(26)
info%f(32) = info%f(26)
info%f(33) = sin1n*sinhn*sinhn/0.0164d0

! Compute V0+u from Schureman (1958)

hpi=pi/2
info%v0_u( 1) = tt-3*s+hp+p+hpi+2*xi-nu	! Q1
info%v0_u( 2) = tt-2*s+hp+hpi+2*xi-nu		! O1
info%v0_u( 3) = tt+hp-hpi-nuprim		! K1
info%v0_u( 4) = 2*(tt-2*s+hp+p+xi-nu)		! 2N2
info%v0_u( 5) = 2*tt-3*s+2*hp+p+2*xi-2*nu	! N2
info%v0_u( 6) = 2*(tt-s+hp+xi-nu)		! M2
info%v0_u( 7) = 2*(tt+hp-nusec)		! K2
info%v0_u( 8) = 2*tt				! S2
info%v0_u( 9) = tt-hp+hpi			! P1
info%v0_u(10) = 4*(tt-s+hp+xi-nu)		! M4
info%v0_u(11) = tt				! S1
info%v0_u(12) = 2*(s-xi)			! Mf
info%v0_u(13) = s-p				! Mm
info%v0_u(14) = 3*s-p-2*xi			! Mtm
info%v0_u(15) = 2*(2*s-hp-xi)			! MSqm
info%v0_u(16) = 2*tt-3*s+4*hp-p+2*xi-2*nu	! Nu2
info%v0_u(17) = 2*tt-4*s+4*hp+2*xi-2*nu	! Mu2
info%v0_u(18) = 2*tt-s+2*hp-p+pi+2*xi-2*nu-r	! L2
info%v0_u(19) = 2*tt-hp+p1			! T2
info%v0_u(20) = 2*tt-5*s+4*hp+p		! Eps2
info%v0_u(21) = 2*tt-s+p+pi+2*xi-2*nu		! Lambda2
info%v0_u(22) = 2*tt+s+2*hp-p-2*nu		! Eta2
info%v0_u(23) = tt-4*s+hp+2*p+hpi+2*xi-nu	! 2Q1
info%v0_u(24) = tt-4*s+3*hp+hpi+2*xi-nu	! Sigma1
info%v0_u(25) = tt-3*s+3*hp-p+hpi+2*xi-nu	! Ro1
info%v0_u(26) = tt-s+hp+p-hpi-nu		! M11
info%v0_u(27) = tt-s+hp-p-hpi+2*xi-nu		! M12
info%v0_u(28) = tt-s+3*hp-p-hpi-nu		! Ki1
info%v0_u(29) = tt-2*hp+p1+hpi			! Pi1
info%v0_u(30) = tt+3*hp-hpi			! Phi1
info%v0_u(31) = tt+s-hp+p-hpi-nu		! Teta1
info%v0_u(32) = tt+s+hp-p-hpi-nu		! J1
info%v0_u(33) = tt+2*s+hp-hpi-2*xi-nu		! OO1

info%v0_u = mod(info%v0_u,2*pi)
end subroutine fes_astron

elemental function expj (x)
real(eightbytereal), intent(in) :: x
complex(eightbytereal) :: expj
expj = dcmplx(cos(x),sin(x))
end function expj

end subroutine festide

!&festideinit -- Initialize FES tide model
!+
subroutine festideinit (name, wantload, info)
character(*), intent(in) :: name
logical, intent(in) :: wantload
type(festideinfo), intent(inout) :: info
!
! Allocate memory for FES tide modeling and read grids into memory.
! When wantload is .TRUE., loading tide grids are loaded and load tide
! will be computed. When .FALSE., festide will return a zero load tide.
!
! Input arguments:
!  name     : Name of the FES tide model (FES95.2.1, FES99, or
!             FES2002, FES2004, or FES2012)
!  wantload : Specify that load tide has to be computed
!
! Output argument:
!  info     : Struct initialised by festideinit
!-
character(256) :: pathname
integer(fourbyteint) ::	istat,lwa

! Determine size and number of grids to load

info%ng = 8
info%nx = 1441
if (name(:5) == 'FES95') then
	info%nx = 721
else if (name(:5) == 'FES99') then
else if (name(:7) == 'FES2002') then
	info%ng = 9
else if (name(:7) == 'FES2004' .or. name(:7) == 'FES2012') then
	info%ng = 9 + 2 + 4
	info%nx = 2881
else
	stop 'Unknown tide model'
endif
info%ny = info%nx/2+1
info%dx = (xmax-xmin)/(info%nx-1)
info%dy = (ymax-ymin)/(info%ny-1)
info%pmin = 0d0

! Allocate memory

if (allocated(info%p_ocean)) stop 'festideinit: Use festidefree first'
allocate (info%p_ocean(info%nx,info%ny,2,info%ng),stat=istat)
if (istat /= 0) stop 'festideinit: not able to allocate memory'

pathname='/user/altim'
call checkenv('ALTIM',pathname,lwa)
pathname(lwa+1:)='/data/'//trim(name)
lwa=len_trim(pathname)

! Ocean tide waves are read from grids.

call fes_read_tide(pathname(:lwa)//'/oceantide.nc',info%p_ocean,info%dz(:,1),info%z0(:,1))

! Read loading grids (when requested)

info%haveload = wantload
if (wantload) then
	allocate (info%p_load(info%nx,info%ny,2,info%ng), stat=istat)
	if (istat /= 0) stop 'festideinit: not able to allocate memory'
	call fes_read_tide(pathname(:lwa)//'/loadtide.nc',info%p_load,info%dz(:,2),info%z0(:,2))
endif

info%t_nodal = 1d30
info%nan = 0d0
info%nan = info%nan/info%nan

contains

!-----------------------------------------------------------------------
! Read amplitude and phase grids from single netCDF file

subroutine fes_read_tide(pathname,work,fact,offs)
use netcdf
character(*), intent(in) :: pathname
real(eightbytereal), intent(out) :: fact(:),offs(:)
integer(twobyteint), intent(out) :: work(:,:,:,:)
integer(fourbyteint) :: ncid,kx,ky,kw
real(eightbytereal) :: dummy(2)
character(16) :: varnm

write (0,600) trim(pathname)
600 format ('(Loading FES tide: ',a,')')

! Open and verify netCDF file

if (nf90_open(pathname,nf90_nowrite,ncid) /= nf90_noerr) stop 'festideinit: Unable to open file'
if (nf90_inquire_dimension(ncid,1,len=kx) /= nf90_noerr) stop 'festideinit: Error in dimensions'
if (nf90_inquire_dimension(ncid,2,len=ky) /= nf90_noerr) stop 'festideinit: Error in dimensions'
if (kx /= info%nx .or. ky /= info%ny) stop 'festideinit: Unexpected grid size'
if (nf90_get_att(ncid,1,'actual_range',dummy) /= nf90_noerr) stop 'festideinit: Error in dimensions'
if (dummy(1) /= xmin .or. dummy(2) /= xmax) stop 'festideinit: Unexpected longitude range'
if (nf90_get_att(ncid,2,'actual_range',dummy) /= nf90_noerr) stop 'festideinit: Error in dimensions'
if (dummy(1) /= ymin .or. dummy(2) /= ymax) stop 'festideinit: Unexpected latitude range'

! Read individual grids

do kw = 1,info%ng
	if (nf90_inquire_variable(ncid,2*kw+1,varnm) /= nf90_noerr) stop 'festideinit: Missing wave'
	if (varnm(5:6) /= wave(kw)(1:2)) stop 'festideinit: Wrong wave'
	if (nf90_get_att(ncid,2*kw+1,'scale_factor',fact(kw)) /= nf90_noerr) fact(kw)=1d0
	if (nf90_get_att(ncid,2*kw+1,'add_offset',offs(kw)) /= nf90_noerr) offs(kw)=0d0
	if (nf90_get_var(ncid,2*kw+1,work(:,:,1,kw)) /= nf90_noerr) stop 'festideinit: Missing wave'
	if (nf90_get_var(ncid,2*kw+2,work(:,:,2,kw)) /= nf90_noerr) stop 'festideinit: Missing wave'
enddo
if (nf90_close(ncid) /= nf90_noerr) return

end subroutine fes_read_tide

end subroutine festideinit

!&festidefree -- Free up space allocated by festideinit
!+
subroutine festidefree (info)
type(festideinfo), intent(inout) :: info
!
! This routine frees up memory allocated by a previous call to festideinit
!
! Argument:
!  info     : Struct initialised by festideinit
!-
if (allocated(info%p_ocean)) deallocate (info%p_ocean)
if (allocated(info%p_load)) deallocate (info%p_load)
end subroutine festidefree

end module tides_fes
