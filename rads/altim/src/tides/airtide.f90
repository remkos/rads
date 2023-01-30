module tides_air
use typesizes
integer(fourbyteint), parameter :: nw=6
real(eightbytereal), parameter :: pi=4d0*atan(1d0),rad=pi/180d0,xmin=0d0,xmax=360d0,ymin=-90d0,ymax=90d0
character(2), parameter :: wave(nw)=(/'S1','P1','K1','S2','T2','R2'/)
type :: airtideinfo
	integer(twobyteint), allocatable :: p_air(:,:,:,:)
	real(eightbytereal) :: dx,dy,z0(nw),dz(nw),nan
	integer(fourbyteint) :: nx,ny
endtype

contains

!*airtide -- Compute air tides according to Ray and Ponte model
!+
function airtide (info, utc, lat, lon)
type(airtideinfo), intent(inout) :: info
real(eightbytereal), intent(in) :: utc, lat, lon
real(eightbytereal) :: airtide
!
! Compute air tide for given time and location based on grids of the
! harmonic constituents S1, P1, K1, S2, T2, and R2 determined by Ray
! and Ponte (2003). Amplitudes and phases for these tides are contained
! in a single netCDF file.
!
! The resulting value TIDE is in millibar and is part of the ocean
! tide model. However, this contribution should be removed from the
! sea level pressure before computing the inverse barometric correction.
!
! The input grids can be found in $ALTIM/data/<tide>, where $ALTIM
! is an environment variable and <tide> is the argument in the
! call to airtide. The grids are in netCDF format and were converted
! using the program 'got2nc'.
!
! To initialize the computation, the function airtideinit should be
! called first. It allocates the appropriate amount of memory and
! loads the grids into memory. To release the memory for further
! use, call airtidefree.
!
! Longitude and latitude are to be specified in degrees; time in UTC
! seconds since 1 Jan 1985. All predicted tides are output in meters.
! If the tide is requested in a point where it is not defined, NaN
! (Not-a-Number) is returned.
!
! Input arguments:
!  info     : Structure initialised by airtideinit
!  utc      : UTC time in seconds since 1 Jan 1985
!  lat      : Latitude (degrees)
!  lon      : Longitude (degrees)
!
! Returned value:
!  airtide  : Air tide (millibars)
!
! Reference:
!  R D Ray and R M Ponte, Barometric tides from ECMWF operational analyses,
!  Annales Geophysicae, 21: 1897-1910, 2003.
!-
! $Log: airtide.f90,v $
! Revision 1.13  2014/02/08 16:39:01  rads
! - Replaced cmplx by dcmplx to avoid loss of precision
!
! Revision 1.12  2011/05/24 14:47:48  rads
! - Created module "tides"
! - Assign "info" argument on initialisation
! - Renamed *free to *tidefree, *init to *tideinit
!
! Revision 1.11  2011/03/08 05:45:28  rads
! - Now read from a single netCDF file airtide.nc
! - Include all six airtide components
!
! Revision 1.8  2008/03/17 17:10:15  rads
! - Ported from Fortran 77 to 90. Speed performance improvement.
! - Fixed scale error by factor 1d-4.
!
! Revision 1.7  2006/09/29 15:08:50  rads
! - Send progress reports to standard error, not standard output
!
! Revision 1.6  2006/09/25 00:45:33  rads
! - Use NetCDF grids as input
!
! Revision 1.4  2006/07/28 17:34:47  rads
! - Avoiding use of %val
!
! Revision 1.3  2005/11/07 22:21:00  rads
! - Use proper time argument
!
! Revision 1.1  2005/10/27 15:20:32  rads
! - Preliminary version of airtide model
!-----------------------------------------------------------------------
real(eightbytereal) :: slon,t,t1,t2,h,pp,arg(nw)
complex(eightbytereal) :: val(nw),w(nw)

! If latitude out of range, bail out

if (lat < ymin .or. lat > ymax) then
	airtide = info%nan
	return
endif

! Limit longitude to interval of grids

slon = lon
if (slon.gt.xmax) then
	slon = slon-360d0
else if (slon.lt.xmin) then
	slon = slon+360d0
endif

! Compute basic astronomical longitudes

t = utc/86400d0
t1 = modulo (t, 1d0) * 360d0
t2 = 2d0 * t1

t = (t - 5478.4993d0) / 36525d0							! Julian centuries since 1.5 Jan 2000
h = mod(280.4662556d0 +  36000.76983081d0 * t, 360d0)	! Mean longitude of sun
pp=     282.94d0      +      1.7192d0     * t			! Mean longitude of solar perigee

! Compute tidal arguments

arg(1) = t1	+ 180d0		! S1 (Doodson's phase)
arg(2) = t1 - h - 90d0	! P1
arg(3) = t1 + h + 90d0	! K1
arg(4) = t2				! S2
arg(5) = t2 - h + pp	! T2
arg(6) = t2 + h - pp	! R2

! Compute all ocean tide components

call air_interp(info%p_air,info%dz,info%z0,slon,lat,val)
w = expj(-arg*rad)
airtide = sum(dble(w*val))

contains

!-----------------------------------------------------------------------
! Interpolate a value from a grid of data at the desired location.
! Interpolation is bilinear.
subroutine air_interp(work,fact,offs,x,y,val)
real(eightbytereal), intent(in) :: x,y,fact(:),offs(:)
complex(eightbytereal), intent(out) :: val(:)
integer(fourbyteint) :: kw,i0,j0,i,j
integer(twobyteint) :: work(:,:,:,:)
real(eightbytereal) :: ptot,xij,yij,pds

val = (0d0,0d0)

xij = (x-xmin)/info%dx+1
yij = (y-ymin)/info%dy+1
i0 = min(info%nx-1,int(xij))
j0 = min(info%ny-1,int(yij))

ptot = 0d0

do i = i0,i0+1
	do j = j0,j0+1
		pds = (1d0-abs(i-xij))*(1d0-abs(j-yij))
		ptot = ptot+pds
		forall (kw=1:nw) val(kw)=val(kw)+(work(i,j,1,kw)*fact(kw)+offs(kw))*expj(work(i,j,2,kw)*1d-2*rad)*pds
	enddo
enddo

val=val/ptot

end subroutine air_interp

elemental function expj (x)
real(eightbytereal), intent(in) :: x
complex(eightbytereal) :: expj
expj = dcmplx(cos(x),sin(x))
end function expj

end function airtide

!&airtideinit -- Initialize air tide model
!+
subroutine airtideinit (name, info)
character(*), intent(in) :: name
type(airtideinfo), intent(inout) :: info
!
! Allocate memory for air tide modeling and read grids into memory.
!
! Input argument:
!  name     : Subdirectory of airtide model in $ALTIM/data
!
! Output argument:
!  info     : Structure initialised by airtideinit
!-
integer(fourbyteint) :: lwa
character(256) :: pathname

! Determine size and number of grids to load

info%nx = 181
info%ny = info%nx/2+1
info%dx = (xmax-xmin)/(info%nx-1)
info%dy = (ymax-ymin)/(info%ny-1)

! Allocate memory

if (allocated(info%p_air)) stop 'airtideinit: Use airtidefree first'
allocate (info%p_air(info%nx,info%ny,2,nw))

pathname = '/user/altim'
call checkenv('ALTIM',pathname,lwa)
pathname(lwa+1:) = '/data/'//trim(name)
lwa=len_trim(pathname)

! Six (nw) waves are read from grid

call air_read_tide(pathname(:lwa)//'/airtide.nc',info%p_air,info%dz,info%z0)

info%nan=0d0
info%nan=info%nan/info%nan

contains

!-----------------------------------------------------------------------
! Read real and imaginary grids in netCDF format

subroutine air_read_tide(pathname,work,fact,offs)
use netcdf
character(*), intent(in) :: pathname
real(eightbytereal), intent(out) :: fact(:),offs(:)
integer(twobyteint), intent(out) :: work(:,:,:,:)
integer(fourbyteint) :: ncid,kx,ky,kw
real(eightbytereal) :: dummy(2)
character(16) :: varnm

write (0,600) trim(pathname)
600 format ('(Loading air tide: ',a,')')

! Open and verify netCDF file

if (nf90_open(pathname,nf90_nowrite,ncid) /= nf90_noerr) stop 'airtideinit: Unable to open file'
if (nf90_inquire_dimension(ncid,1,len=kx) /= nf90_noerr) stop 'airtideinit: Error in dimensions'
if (nf90_inquire_dimension(ncid,2,len=ky) /= nf90_noerr) stop 'airtideinit: Error in dimensions'
if (kx /= info%nx .or. ky /= info%ny) stop 'airtideinit: Unexpected grid size'
if (nf90_get_att(ncid,1,'actual_range',dummy) /= nf90_noerr) stop 'airtideinit: Error in dimensions'
if (dummy(1) /= xmin .or. dummy(2) /= xmax) stop 'airtideinit: Unexpected longitude range'
if (nf90_get_att(ncid,2,'actual_range',dummy) /= nf90_noerr) stop 'airtideinit: Error in dimensions'
if (dummy(1) /= ymin .or. dummy(2) /= ymax) stop 'airtideinit: Unexpected latitude range'

! Read individual grids

do kw = 1,nw
	if (nf90_inquire_variable(ncid,2*kw+1,varnm) /= nf90_noerr) stop 'airtideinit: Missing wave'
	if (varnm(5:6) /= wave(kw)(1:2)) stop 'airtideinit: Wrong wave'
	if (nf90_get_att(ncid,2*kw+1,'scale_factor',fact(kw)) /= nf90_noerr) fact(kw)=1d0
	if (nf90_get_att(ncid,2*kw+1,'add_offset',offs(kw)) /= nf90_noerr) offs(kw)=0d0
	if (nf90_get_var(ncid,2*kw+1,work(:,:,1,kw)) /= nf90_noerr) stop 'airtideinit: Missing wave'
	if (nf90_get_var(ncid,2*kw+2,work(:,:,2,kw)) /= nf90_noerr) stop 'airtideinit: Missing wave'
enddo
if (nf90_close(ncid) /= nf90_noerr) return

end subroutine air_read_tide

end subroutine airtideinit

!&airtidefree -- Free up space allocated by GOTINIT
!+
subroutine airtidefree (info)
type(airtideinfo), intent(inout) :: info
!
! This routine frees up memory allocated by a previous call to GOTINIT
!
! Argument:
!  info     : Structure initialised by airtideinit
!-
if (allocated(info%p_air)) deallocate (info%p_air)
end subroutine airtidefree

end module tides_air
