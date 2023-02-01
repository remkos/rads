module tides_got
use typesizes
integer(twobyteint), parameter, private :: undef=-32768
integer(fourbyteint), parameter, private :: nw=28
real(eightbytereal), parameter, private :: pi=4d0*atan(1d0),rad=pi/180d0,xmin=0d0,xmax=360d0,ymin=-90d0,ymax=90d0
character(8), parameter, private :: wave(nw) = (/ &
	'Q1      ','O1      ','P1      ','K1      ','N2      ','M2      ','S2      ', &
	'K2      ','M4      ','S1      ','2Q1     ','sigma1  ','rho1    ','M1      ', &
	'M1      ','chi1    ','pi1     ','phi1    ','theta1  ','J1      ','OO1     ', &
	'2N2     ','mu2     ','nu2     ','lambda2 ','L2      ','L2      ','T2      ' /)
type :: gottideinfo
	integer(fourbyteint) :: nx,ny,ng
	logical :: haveload
	integer(twobyteint), allocatable :: p_ocean(:,:,:,:),p_load(:,:,:,:)
	real(eightbytereal) :: dx,dy,z0(nw,2),dz(nw,2),f(nw),u(nw),t_nodal,nan,pmin
endtype

contains

!*gottide -- Compute tides according to GOT model
!+
subroutine gottide (info, utc, lat, lon, tide, tide_load)
use tides_aux
type(gottideinfo), intent(inout) :: info
real(eightbytereal), intent(in) :: utc, lat, lon
real(eightbytereal), intent(out) :: tide, tide_load
!
! Compute ocean tidal height for given time and location from grids
! of harmonic coefficients of one of the GOT models.
! Current version uses the 8 largest constituents in the
! semidiurnal & diurnal bands, with other tides inferred.
! This routine is heavily based on the routine PERTH3 by Richard Ray.
!
! The input files can be found in $ALTIM/data/<name>, where $ALTIM
! is an environment variable and <name> is the argument in the
! call to gottideinit. The grids are in NetCDF format and were converted
! using the program 'got2nc'.
!
! To initialize the computation, the function gottideinit should be
! called first. It allocates the appropriate amount of memory and
! loads the grids into memory. To release the memory for further
! use, call gottidefree.
!
! Longitude and latitude are to be specified in degrees; time in UTC
! seconds since 1 Jan 1985. All predicted tides are output in meters.
! If the tide is requested in a point where it is not defined, NaN
! (Not-a-Number) is returned.
!
! Note that this routine does NOT compute the long-period equilibrium
! tides. Use the routine lpetide to compute those.
!
! The computation of nodal periods and mean longitudes is optimized
! for the period 1980-2020.
!
! Input arguments:
!  info     : Structure initialised by gottideinit
!  utc      : UTC time in seconds since 1 Jan 1985
!  lat      : Latitude (degrees)
!  lon      : Longitude (degrees)
!
! Output arguments:
!  tide     : Predicted short-period tide (m)
!  tide_load: Predicted loading effect (m)
!-
! $Log: gottide.f90,v $
! Revision 1.32  2019/03/08 08:35:04  rads
! - Return NaN when time stamp is NaN
!
! Revision 1.31  2017/10/12 09:17:55  rads
! - Add generic mean_longtitudes routine
!
! Revision 1.30  2014/05/03 16:21:46  rads
! - Allow changing of interpolation limit using info%pmin
!
! Revision 1.29  2014/02/08 16:39:01  rads
! - Replaced cmplx by dcmplx to avoid loss of precision
!
! Revision 1.28  2011/08/12 20:38:15  rads
! - Updated astronomical arguments u and f every second, instead of every day
! - No need to store info%arg, just arg will do (better parallelisation)
!
! Revision 1.27  2011/05/24 14:47:49  rads
! - Created module "tides"
! - Assign "info" argument on initialisation
! - Renamed *free to *tidefree, *init to *tideinit
!
! Revision 1.26  2009/05/20 15:53:00  rads
! - Replaced sum(a*b) by dot_product(a,b) (no impact on results)
!
! Revision 1.25  2009/05/12 20:55:39  rads
! - All tidal components now stored in one grid
!
! Revision 1.24  2009/01/16 16:35:34  rads
! - Arrays of characters need to be initialized with same length (gfortran 4.4)
!
! Revision 1.23  2008/06/03 21:30:43  rads
! - No longer computes long-period tides
!
! Revision 1.22  2008/03/17 23:38:43  rads
! - Marginal adjustments in mean longitudes to center on year 2000.
! - Adapted to GOT4.x (includes S1 and M4 tides) based on perth3.f by Richard Ray.
!
! Revision 1.21  2008/03/17 17:20:21  rads
! - Ported from Fortran 77 to 90. Speed performace improvements.
!
! Revision 1.20  2006/09/29 15:08:50  rads
! - Send progress reports to standard error, not standard output
!
! Revision 1.19  2006/09/25 00:14:40  rads
! - Use netCDF input files
!
! Revision 1.17  2006/08/07 17:24:59  rads
! - Removing obsolete Fortran code so it can be used with gfortran
!
! Revision 1.14  2005/10/27 20:15:16  rads
! - Polished code; no effect on results
!
! Revision 1.12  2005/05/31 20:52:45  rads
! - Do all internal computation in double precision
!
! Revision 1.10  2004/11/23 01:55:03  remko
! - Automatically allocates memory using gottideinit and gottidefree
! - Change of arguments in gottide
!
! 10-Jun-2003 - Allow grids with various scales; increased update of
!               astronomicals to once per day (from once per month)
!  9-Oct-2001 - Adapted by Remko Scharroo for DEOS
! 13-Jan-1998 - Version of PERTH2 by Richard Ray
!-----------------------------------------------------------------------
real(eightbytereal) :: slon,arg(nw)
real(eightbytereal), parameter :: mjd85=46066d0
integer(fourbyteint) :: istat1,istat2
complex(eightbytereal) :: val(nw),w(nw)

! When called with invalid time, reset time reference and return

if (.not.(utc < 1d20)) then
	info%t_nodal=1d30
	tide=info%nan
	tide_load=info%nan
	return
endif

! If latitude out of range, bail out

if (lat < ymin .or. lat > ymax) then
	tide=info%nan
	tide_load=info%nan
	return
endif

! Limit longitude to interval of grids

slon=lon
if (slon > xmax) then
	slon=slon-360d0
else if (slon < xmin) then
	slon=slon+360d0
endif

! Compute basic astronomical mean longitudes

call got_astron(utc)

! Compute all ocean tide components

call got_interp(info%p_ocean,info%dz(:,1),info%z0(:,1),slon,lat,val,istat1)
if (istat1 == 0) then
	tide=info%nan
else
	w = expj(-(arg+info%u)*rad)
	tide=dot_product(info%f,dble(w*val))
endif

! Compute all load tide components

tide_load=0d0
if (.not.info%haveload) return

call got_interp(info%p_load,info%dz(:,2),info%z0(:,2),slon,lat,val,istat2)
if (istat2 == 0) then
	tide_load=info%nan
else
	if (istat1 == 0) w = expj(-(arg+info%u)*rad)
	tide_load=dot_product(info%f,dble(w*val))
endif

contains

!-----------------------------------------------------------------------
! Compute the basic astronomical mean longitudes  s, h, p, N.

subroutine got_astron (time)
real(eightbytereal), intent(in) :: time
!
! These formulae are for the period 1980 - 2020, and were derived
! from the ASTRO5 code in perth3.f by Richard Ray and are based on:
! Jean Meeus, Astronomical Algorithms, 2nd ed., 1998.
!
! Note 1: This routine uses time in UT and does not distinguish
! between the subtle differences of UTC, UT1, etc. However, this is
! more than accurate enough for our purposes.
! The formalae for the mean longitudes depend on dynamic time (DT).
! This routine assumes DT - UT = 63.48 seconds (2000).
!
! Note 2: The argument w is -n', i.e. w is decreasing with time.
!
! TIME is UTC in seconds since 1985.
!
real(eightbytereal) :: t1,t2,s,h,p,w,pp,sinn,sin2n,cosn,cos2n,t

t = time/86400d0
t1 = modulo (t, 1d0) * 360d0
t2 = 2d0*t1

call mean_longitudes (t, s, h, p, w, pp)
w = -w ! Change n' into w = -n'

arg( 1) = t1 + h - 3d0*s + p - 90d0	! Q1
arg( 2) = t1 + h - 2d0*s - 90d0		! O1
arg( 3) = t1 - h - 90d0			! P1
arg( 4) = t1 + h + 90d0			! K1
arg( 5) = t2 + 2d0*h - 3d0*s + p	! N2
arg( 6) = t2 + 2d0*h - 2d0*s		! M2
arg( 7) = t2				! S2
arg( 8) = t2 + 2d0*h			! K2
arg( 9) = 2d0*arg(6)			! M4
arg(10) = t1 + 180d0			! S1 (Doodson's phase)
arg(11) = t1 - 4d0*s + h + 2d0*p - 90d0	! 2Q1
arg(12) = t1 - 4d0*s + 3d0*h - 90d0	! sigma1
arg(13) = t1 - 3d0*s + 3d0*h - p - 90d0	! rho1
arg(14) = t1 - s + h - p + 90d0		! M1
arg(15) = t1 - s + h + p + 90d0		! M1
arg(16) = t1 - s + 3d0*h - p + 90d0	! chi1
arg(17) = t1 - 2d0*h + pp - 90d0	! pi1
arg(18) = t1 + 3d0*h + 90d0		! phi1
arg(19) = t1 + s - h + p + 90d0		! theta1
arg(20) = t1 + s + h - p + 90d0		! J1
arg(21) = t1 + 2d0*s + h + 90d0		! OO1
arg(22) = t2 - 4d0*s + 2d0*h + 2d0*p	! 2N2
arg(23) = t2 - 4d0*s + 4d0*h		! mu2
arg(24) = t2 - 3d0*s + 4d0*h - p	! nu2
arg(25) = t2 - s + p + 180d0		! lambda2
arg(26) = t2 - s + 2d0*h - p + 180d0	! L2
arg(27) = t2 - s + 2d0*h + p		! L2
arg(28) = t2 - h + pp			! T2

! Determine nodal corrections f and u (only when more than 1 day is passed since last time)

if (abs(time-info%t_nodal) > 86400d0) then
	info%t_nodal = time
	sinn = sin(w*rad)
	cosn = cos(w*rad)
	sin2n = sin(2*w*rad)
	cos2n = cos(2*w*rad)

	info%f( 1) = 1.009d0 + 0.187d0*cosn - 0.015d0*cos2n
	info%f( 2) = info%f(1)
	info%f( 4) = 1.006d0 + 0.115d0*cosn - 0.009d0*cos2n
	info%f( 5) = 1.000d0 - 0.037d0*cosn
	info%f( 6) = info%f(5)
	info%f( 8) = 1.024d0 + 0.286d0*cosn + 0.008d0*cos2n
	info%f( 9) = info%f(6)*info%f(6)
	info%f(11) = sqrt((1.0d0 + 0.189d0*cosn - 0.0058d0*cos2n)**2 + (0.189d0*sinn - 0.0058d0*sin2n)**2)
	info%f(12) = info%f(11)
	info%f(13) = info%f(11)
	info%f(14) = sqrt((1.0d0 + 0.185d0*cosn)**2 + (0.185d0*sinn)**2)
	info%f(15) = sqrt((1.0d0 + 0.201d0*cosn)**2 + (0.201d0*sinn)**2)
	info%f(16) = sqrt((1.0d0 + 0.221d0*cosn)**2 + (0.221d0*sinn)**2)
	info%f(20) = sqrt((1.0d0 + 0.198d0*cosn)**2 + (0.198d0*sinn)**2)
	info%f(21) = sqrt((1.0d0 + 0.640d0*cosn + 0.134d0*cos2n)**2 + (0.640d0*sinn + 0.134d0*sin2n)**2)
	info%f(22) = sqrt((1.0d0 - 0.0373d0*cosn)**2 + (0.0373d0*sinn)**2)
	info%f(23) = info%f(22)
	info%f(24) = info%f(22)
	info%f(26) = info%f(22)
	info%f(27) = sqrt((1.0d0 + 0.441d0*cosn)**2 + (0.441d0*sinn)**2)

	info%u( 1) = 10.8d0*sinn - 1.3d0*sin2n
	info%u( 2) = info%u(1)
	info%u( 4) = -8.9d0*sinn + 0.7d0*sin2n
	info%u( 5) = -2.1d0*sinn
	info%u( 6) = info%u(5)
	info%u( 8) = -17.7d0*sinn + 0.7d0*sin2n
	info%u( 9) = 2d0*info%u(6)
	info%u(11) = atan2(0.189d0*sinn - 0.0058d0*sin2n, 1.0d0 + 0.189d0*cosn - 0.0058d0*sin2n)/rad
	info%u(12) = info%u(11)
	info%u(13) = info%u(11)
	info%u(14) = atan2( 0.185d0*sinn, 1.0d0 + 0.185d0*cosn)/rad
	info%u(15) = atan2(-0.201d0*sinn, 1.0d0 + 0.201d0*cosn)/rad
	info%u(16) = atan2(-0.221d0*sinn, 1.0d0 + 0.221d0*cosn)/rad
	info%u(20) = atan2(-0.198d0*sinn, 1.0d0 + 0.198d0*cosn)/rad
	info%u(21) = atan2(-0.640d0*sinn - 0.134d0*sin2n, 1.0d0 + 0.640d0*cosn + 0.134d0*cos2n)/rad
	info%u(22) = atan2(-0.0373d0*sinn, 1.0d0 - 0.0373d0*cosn)/rad
	info%u(23) = info%u(22)
	info%u(24) = info%u(22)
	info%u(26) = info%u(22)
	info%u(27) = atan2(-0.441d0*sinn, 1.0d0 + 0.441d0*cosn)/rad
endif
end subroutine got_astron

!-----------------------------------------------------------------------
! Interpolates a value from a grid of data at the desired location.
! Interpolation is bilinear.

subroutine got_interp(work,fact,offs,x,y,val,istat)
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
! Check if the major constituents are available:
		do kw = 1,info%ng
			if (work(i,j,1,kw) == undef) cycle jloop
		enddo

		istat = istat+1
		pds = (1d0-abs(i-xij))*(1d0-abs(j-yij))
		ptot = ptot+pds

		forall (kw=1:info%ng) val(kw)=val(kw)+(work(i,j,1,kw)*fact(kw)+offs(kw))*expj(work(i,j,2,kw)*1d-2*rad)*pds
	enddo jloop
enddo

if (istat == 0 .or. ptot <= info%pmin) then
	istat = 0
	return
endif
val(1:info%ng) = val(1:info%ng)/ptot

! Infer additional constituents by admittance

val(11) =  0.2630d0*val(1) - 0.0252d0*val(2)	! 2Q1
val(12) =  0.2970d0*val(1) - 0.0264d0*val(2)	! sigma1
val(13) =  0.1640d0*val(1) + 0.0048d0*val(2)	! rho1
val(14) =  0.0140d0*val(2) + 0.0101d0*val(4)	! M1
val(15) =  0.0389d0*val(2) + 0.0282d0*val(4)	! M1
val(16) =  0.0064d0*val(2) + 0.0060d0*val(4)	! chi1
val(17) =  0.0030d0*val(2) + 0.0171d0*val(4)	! pi1
val(18) = -0.0015d0*val(2) + 0.0152d0*val(4)	! phi1
val(19) = -0.0065d0*val(2) + 0.0155d0*val(4)	! theta1
val(20) = -0.0389d0*val(2) + 0.0836d0*val(4)	! J1
val(21) = -0.0431d0*val(2) + 0.0613d0*val(4)	! OO1
val(22) =  0.2640d0*val(5) - 0.0253d0*val(6)	! 2N2
val(23) =  0.2980d0*val(5) - 0.0264d0*val(6)	! mu2
val(24) =  0.1650d0*val(5) + 0.00487d0*val(6)	! nu2
val(25) =  0.0040d0*val(6) + 0.0074d0*val(7)	! lambda2
val(26) =  0.0131d0*val(6) + 0.0326d0*val(7)	! L2
val(27) =  0.0033d0*val(6) + 0.0082d0*val(7)	! L2
val(28) =  0.0585d0*val(7)						! T2
end subroutine got_interp

elemental function expj (x)
real(eightbytereal), intent(in) :: x
complex(eightbytereal) :: expj
expj = dcmplx(cos(x),sin(x))
end function expj

end subroutine gottide

!&gottideinit -- Initialize GOT tide model
!+
subroutine gottideinit (name, wantload, info)
character(*), intent(in) :: name
logical, intent(in) :: wantload
type(gottideinfo), intent(inout) :: info
!
! Allocate memory for GOT tide modeling and read grids into memory.
! When wantload is .TRUE., loading tide grids are loaded and load tide
! will be computed. When .FALSE., gottide will return a zero load tide.
!
! Input arguments:
!  name     : Name of the GOT tide model (GOT99.2b or GOT00.2)
!  wantload : Specify that load tide has to be computed
!
! Output argument:
!  info     : Structure initialised by gottideinit
!-
integer(fourbyteint) :: istat,lwa
character(256) :: pathname

! Determine size and number of grids to load

info%ng = 8
if (name(:4) == 'GOT4') info%ng = 10
info%nx = 721
info%ny = info%nx/2+1
info%dx = (xmax-xmin)/(info%nx-1)
info%dy = (ymax-ymin)/(info%ny-1)
info%pmin = 0.5d0

! Allocate memory

if (allocated(info%p_ocean)) stop 'gottideinit: Use gottidefree first'
allocate (info%p_ocean(info%nx,info%ny,2,info%ng), stat=istat)
if (istat /= 0) stop 'gottideinit: not able to allocate memory'

pathname='/user/altim'
call checkenv('ALTIM',pathname,lwa)
pathname(lwa+1:) = '/data/'//trim(name)
lwa = len_trim(pathname)

! Major (info%ng) waves are read from grids and 18 are inferred by admittance.

call got_read_tide(pathname(:lwa)//'/oceantide.nc',info%p_ocean,info%dz(:,1),info%z0(:,1))

! Read loading grids (when requested)

info%haveload = wantload
if (wantload) then
	allocate (info%p_load(info%nx,info%ny,2,info%ng), stat=istat)
	if (istat /= 0) stop 'gottideinit: not able to allocate memory'
	call got_read_tide(pathname(:lwa)//'/loadtide.nc',info%p_load,info%dz(:,2),info%z0(:,2))
endif

info%f = 1d0
info%u = 0d0
info%t_nodal = 1d30
info%nan = 0d0
info%nan = info%nan / info%nan

contains

!-----------------------------------------------------------------------
! Read amplitude and phase grids from single netCDF file

subroutine got_read_tide(pathname,work,fact,offs)
use netcdf
character(*), intent(in) :: pathname
real(eightbytereal), intent(out) :: fact(:),offs(:)
integer(twobyteint), intent(out) :: work(:,:,:,:)
integer(fourbyteint) :: ncid,kx,ky,kw
real(eightbytereal) :: dummy(2)
character(16) :: varnm

write (0,600) trim(pathname)
600 format ('(Loading GOT tide: ',a,')')

! Open and verify netCDF file

if (nf90_open(pathname,nf90_nowrite,ncid) /= nf90_noerr) stop 'gottideinit: Unable to open file'
if (nf90_inquire_dimension(ncid,1,len=kx) /= nf90_noerr) stop 'gottideinit: Error in dimensions'
if (nf90_inquire_dimension(ncid,2,len=ky) /= nf90_noerr) stop 'gottideinit: Error in dimensions'
if (kx /= info%nx .or. ky /= info%ny) stop 'gottideinit: Unexpected grid size'
if (nf90_get_att(ncid,1,'actual_range',dummy) /= nf90_noerr) stop 'gottideinit: Error in dimensions'
if (dummy(1) /= xmin .or. dummy(2) /= xmax) stop 'gottideinit: Unexpected longitude range'
if (nf90_get_att(ncid,2,'actual_range',dummy) /= nf90_noerr) stop 'gottideinit: Error in dimensions'
if (dummy(1) /= ymin .or. dummy(2) /= ymax) stop 'gottideinit: Unexpected latitude range'

! Read individual grids

do kw = 1,info%ng
	if (nf90_inquire_variable(ncid,2*kw+1,varnm) /= nf90_noerr) stop 'gottideinit: Missing wave'
	if (varnm(5:6) /= wave(kw)(1:2)) stop 'gottideinit: Wrong wave'
	if (nf90_get_att(ncid,2*kw+1,'scale_factor',fact(kw)) /= nf90_noerr) fact(kw)=1d0
	if (nf90_get_att(ncid,2*kw+1,'add_offset',offs(kw)) /= nf90_noerr) offs(kw)=0d0
	if (nf90_get_var(ncid,2*kw+1,work(:,:,1,kw)) /= nf90_noerr) stop 'gottideinit: Missing wave'
	if (nf90_get_var(ncid,2*kw+2,work(:,:,2,kw)) /= nf90_noerr) stop 'gottideinit: Missing wave'
enddo
if (nf90_close(ncid) /= nf90_noerr) return

end subroutine got_read_tide

end subroutine gottideinit

!&gottidefree -- Free up space allocated by gottideinit
!+
subroutine gottidefree (info)
type(gottideinfo), intent(inout) :: info
!
! This routine frees up memory allocated by a previous call to gottideinit
!
! Argument:
!  info : Struct initialised by gottideinit
!-
if (allocated(info%p_ocean)) deallocate (info%p_ocean)
if (allocated(info%p_load)) deallocate (info%p_load)
end subroutine gottidefree

end module tides_got
