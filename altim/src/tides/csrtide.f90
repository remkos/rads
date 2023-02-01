module tides_csr
use typesizes
save
integer(fourbyteint), parameter :: nx=721,ny=361,ng=12,undef=-2147483647-1
real(eightbytereal), parameter :: pi=4d0*atan(1d0),rad=pi/180d0,xmin=0d0,xmax=360d0,ymin=-90d0,ymax=90d0, &
	dx=0.5d0,dy=0.5d0,fact=1d-3
real(eightbytereal), parameter :: phc(4)=(/290.21d0, 280.12d0, 274.35d0, 343.51d0/), &
	dpd(4)=(/13.1763965d0,0.9856473d0,0.1114041d0,0.0529539d0/), &
	u00(2)=(/0.0298d0,0.0200d0/),u20(2)=(/0.1408d0,0.0905d0/), &
	u21(2)=(/0.0805d0,0.0638d0/),u40(2)=(/0.6002d0,0.3476d0/), &
	u41(2)=(/0.3025d0,0.1645d0/),v41(2)=(/0.1517d0,0.0923d0/), &
	tc=63072d3	! tc = 1.0 Jan 1987
integer(fourbyteint), parameter :: indx(30,2,4)=reshape((/ &	!(  l  ,m,n)
  -3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, &	!( 1:15,1,1)
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, &	!(16:30,1,1)
  -3,-3,-2,-2,-2,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, &	!( 1:15,2,1)
   0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, &	!(16:30,2,1)
   0, 2, 0, 0, 2, 2, 0, 0, 0, 0, 2, 0, 0, 0, 2, &	!( 1:15,1,2)
  -3,-2,-2, 0, 0, 0, 1, 2,-2, 0, 0,-2, 0, 0, 0, &	!(16:30,1,2)
   0, 2, 0, 2, 3,-1, 0, 0, 1, 2, 3,-2,-1, 0, 0, &	!( 1:15,2,2)
   1,-2, 0, 0, 0,-3,-2,-1, 0, 0, 0, 0,-2, 0, 0, &	!(16:30,2,2)
   2, 0, 1, 1,-1,-1, 0, 0, 0, 2, 0,-1, 1, 1,-1, &	!( 1:15,1,3)
   0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0,-2, 0, 0, &	!(16:30,1,3)
   3, 1, 2, 0, 0, 1, 1, 1, 1,-1,-1, 2, 0, 0, 0, &	!( 1:15,2,3)
   0, 1,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, &	!(16:30,2,3)
   0, 0,-1, 0,-1, 0,-2,-1, 0, 0, 0, 0, 0, 1, 0, &	!( 1:15,1,4)
   0,-1, 0,-1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, &	!(16:30,1,4)
   0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0,-1, 0, &	!( 1:15,2,4)
   0, 0, 0, 0, 1, 0, 0, 0,-1, 0, 1, 2, 0, 0, 1  &	!(16:30,2,4)
   /),(/30,2,4/))
real(eightbytereal), parameter :: pha(30,2)=reshape((/ &				!(  l  ,m)
   90d0, 90d0, 90d0, 90d0, 90d0, 90d0,270d0, 90d0, 90d0,270d0, &	!( 1:10,1)
  270d0,270d0,270d0,270d0,270d0, 13d0,270d0, 90d0, 90d0,270d0, &	!(11:20,1)
  270d0,-13d0,270d0,270d0,270d0,270d0,270d0,270d0,270d0,270d0, &	!(21:30,1)
    0d0,  0d0,  0d0,  0d0, 77d0,103d0,180d0,  0d0, 77d0,  0d0, &	!( 1:10,2)
   77d0,180d0,103d0,180d0,  0d0, 77d0,180d0,180d0,  0d0,  0d0, &	!(11:20,2)
  283d0,  0d0,262d0,180d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0  &	!(21:30,2)
  /),(/30,2/))
real(eightbytereal), parameter :: amp(30,2)=reshape((/ &				!(  l  ,m)
  0.00663d0,0.00802d0,0.00947d0,0.05019d0,0.00180d0,0.00954d0, &	!( 1: 6,1)
  0.00152d0,0.04946d0,0.26218d0,0.00171d0,0.00343d0,0.00741d0, &	!( 7:12,1)
  0.02062d0,0.00414d0,0.00394d0,0.00713d0,0.00137d0,0.12200d0, &	!(13:18,1)
  0.00730d0,0.36874d0,0.05002d0,0.00293d0,0.00524d0,0.00395d0, &	!(19:24,1)
  0.02061d0,0.00409d0,0.00342d0,0.00169d0,0.01128d0,0.00723d0, &	!(25:30,1)
  0.00180d0,0.00467d0,0.01601d0,0.01932d0,0.00130d0,0.00102d0, &	!( 1: 6,2)
  0.00451d0,0.12100d0,0.00113d0,0.02298d0,0.00106d0,0.00190d0, &	!( 7:12,2)
  0.00218d0,0.02358d0,0.63194d0,0.00193d0,0.00466d0,0.01786d0, &	!(13:18,2)
  0.00447d0,0.00197d0,0.01718d0,0.29402d0,0.00305d0,0.00102d0, &	!(19:24,2)
  0.07994d0,0.02382d0,0.00259d0,0.00086d0,0.00447d0,0.00195d0  &	!(25:30,2)
  /),(/30,2/))

type :: csrtideinfo
	integer(fourbyteint), allocatable :: p_ocean(:,:,:),p_load(:,:,:)
	real(eightbytereal) :: nan, t_nodal
	logical :: pseudo,radiation,haveload
	real(eightbytereal) :: ct(30,2),st(30,2)
endtype

contains

!*csrtide -- Compute tides according to CSR model
!+
subroutine csrtide (info, utc, lat, lon, tide, tide_load)
type(csrtideinfo), intent(inout) :: info
real(eightbytereal), intent(in) :: utc, lat, lon
real(eightbytereal), intent(out) :: tide, tide_load
!
! Compute ocean tidal height for given time and location from grids
! of ortho weights for one of the CSR models.
! This routine is heavily based on the routine CSRTPTIDE by Richard Eanes.
!
! For partical purposes I have changed the input to netCDF grids.
! These grids can be found in $ALTIM/data/csr_tides, where $ALTIM is an
! environment variable.
!
! To initialize the computation, the function csrtideinit should be
! called first. It allocates the appropriate amount of memory and
! loads the grids into memory. To release the memory for further
! use, call csrtidefree.
!
! Longitude and latitude are to be specified in degrees; time in UTC
! seconds since 1 Jan 1985. All predicted tides are output in meters.
! If the tide is requested in a point where it is not defined, NaN
! (Not-a-Number) is returned.
!
! Note that the tidal prediction provided by this routine does NOT
! include long-period tides. Use the routine lpetide to compute those.
!
! Input arguments:
!  info     : Structure initialised by csrtideinit
!  utc      : UTC time in seconds since 1 Jan 1985
!  lat      : Latitude (degrees)
!  lon      : Longitude (degrees)
!
! Output arguments:
!  tide     : Predicted short-period tide (m)
!  tide_load: Predicted loading effect (m)
!-
! $Log: csrtide.f90,v $
! Revision 1.11  2014/02/08 16:39:01  rads
! - Replaced cmplx by dcmplx to avoid loss of precision
!
! Revision 1.10  2011/05/24 14:47:48  rads
! - Created module "tides"
! - Assign "info" argument on initialisation
! - Renamed *free to *tidefree, *init to *tideinit
!
! Revision 1.9  2008/06/03 21:28:00  rads
! - No longer computes long-period equilibrium tide
!
! Revision 1.8  2008/03/17 17:25:21  rads
! - Ported from Fortran 77 to 90
!
! Revision 1.6  2006/07/28 17:34:47  rads
! - Avoiding use of %val
!
! Revision 1.5  2006/06/19 13:01:21  rads
! - Removed superfluous variables; streamlined code further
!
! Revision 1.4  2006/06/07 21:31:06  rads
! - Final completely remodelled routine using NetCDF data files
!-----------------------------------------------------------------------
real(eightbytereal) :: ruv(ng),slon

! When called with invalid time, simply ignore

if (utc > 1d20) return

! If latitude out of range, bail out

if (lat < ymin .or. lat > ymax) then
	tide = info%nan
	tide_load = info%nan
	return
endif

! Limit longitude to interval of grids

slon=lon
if (slon > xmax+dx) then
	slon = slon-360d0
else if (slon < xmin) then
	slon = slon+360d0
endif

! Call CSR tide model routines

if (csr_interp(info%p_ocean,slon,lat,ruv)) then
	call csrweights(ruv,utc,tide)
else
	tide = info%nan
endif
if (info%haveload .and. csr_interp(info%p_load,slon,lat,ruv)) then
	call csrweights(ruv,utc,tide_load)
else
	tide_load = 0d0
endif

contains

function csr_interp(uv,x,y,ruv)
integer(fourbyteint), intent(in) :: uv(:,:,:)
real(eightbytereal), intent(in) :: x,y
real(eightbytereal), intent(out) :: ruv(:)
real(eightbytereal) :: pds,xij,yij,ptot
integer(fourbyteint) :: i,j,i0,j0,kw
logical :: csr_interp

! Compute indices for desired position

xij = (x-xmin)/dx+1
yij = (y-ymin)/dy+1
i0 = min(nx-1,int(xij))
j0 = min(ny-1,int(yij))

! Sum up the weighted values for the corners

ptot = 0d0
ruv = 0d0

do i=i0,i0+1
	do j=j0,j0+1
		pds = (1d0-abs(i-xij))*(1d0-abs(j-yij))
		ptot = ptot+pds
		if (uv(i,j,1) /= undef) forall (kw=1:ng) ruv(kw)=ruv(kw)+uv(i,j,kw)*pds
	enddo
enddo

if (ptot < 0.2d0) then
	csr_interp=.false.
else
	ruv = ruv/ptot*fact
	csr_interp=.true.
endif
end function csr_interp

! Compute the tide at time ts from the orthoweights u,v.
! This is a kernel routine, meant to be called only by tptide.

subroutine csrweights(uv,ts,tide)
real(eightbytereal), intent(in) :: uv(2,3,2),ts
real(eightbytereal), intent(out) :: tide
real(eightbytereal) :: shpn(4),phi,td,e
complex(eightbytereal) :: ab(3),pq(3,2),rp
integer(fourbyteint) ::	l,m

if (ts /= info%t_nodal) then
	td = (ts - tc)/86400d0
! Compute 4 principal mean longitudes for a given time td
	shpn = phc + td*dpd
	e = 360d0 * mod(td,1d0) - shpn(1) + shpn(2)
	do m=1,2
! Compute complex orthotides pq(3,2) for given time ts in terms
! of potential complex amplitudes ab(3)
	ab = (0d0,0d0)
		do l=1,30
			phi = m*e + pha(l,m) + sum(indx(l,m,:)*shpn(:))
			rp = amp(l,m)*expj(-phi*rad)
			ab(1) = ab(1) + rp
			ab(2) = ab(2) + rp*info%ct(l,m)
			ab(3) = ab(3) + rp*info%st(l,m)
		enddo
		pq(1,m) = u00(m)*ab(1)
		pq(2,m) = u20(m)*ab(1) - u21(m)*ab(2)
		pq(3,m) = u40(m)*ab(1) - u41(m)*ab(2) + v41(m)*ab(3)
	enddo
endif

! Compute tide

tide = 0d0
do m=1,2
	tide = tide + sum(uv(1,:,m)*dble(pq(:,m)) + uv(2,:,m)*dimag(pq(:,m)))
enddo

info%t_nodal = ts

end subroutine csrweights

elemental function expj (x)
real(eightbytereal), intent(in) :: x
complex(eightbytereal) :: expj
expj = dcmplx(cos(x),sin(x))
end function expj

end subroutine csrtide

!&csrtideinit -- Initialize CSR tide model
!+
subroutine csrtideinit (name, wantload, info)
character(*), intent(in) :: name
logical, intent(in) :: wantload
type(csrtideinfo), intent(inout) :: info
!
! Allocate memory for CSR tide modeling and read grids into memory.
! When WANTLOAD is .TRUE., loading tide grids are loaded and load tide
! will be computed. When .FALSE., csrtide will return a zero load tide.
!
! Input arguments:
!  name     : Name of the CSR tide model (csr_3.0 or csr_4.0)
!  wantload : Specify that load tide has to be computed
!
! Output argument:
!  info     : Structure initialised by csrtideinit
!-
integer(fourbyteint) :: lwa,istat,l,m
character(160) :: pathname
real(eightbytereal) :: theta

if (allocated(info%p_ocean).or.allocated(info%p_load)) stop "csrtideinit: Use csrtidefree first"
allocate (info%p_ocean(nx,ny,ng),stat=istat)
if (istat /= 0) stop "csrtideinit: not able to allocate memory"
if (wantload) allocate (info%p_load(nx,ny,ng),stat=istat)
if (istat /= 0) stop "csrtideinit: not able to allocate memory"

info%haveload = wantload
pathname = '/user/altim'
call checkenv('ALTIM',pathname,lwa)
pathname(lwa+1:) = '/data/csr_tides/ortho_'//trim(name)
lwa = len_trim(pathname)
if (.not.csr_read_tide(pathname(:lwa)//".ot.nc",info%p_ocean)) stop "csrtideinit: Error opening file."
if (info%haveload) info%haveload=csr_read_tide(pathname(:lwa)//".rload.nc",info%p_load)

info%nan = 0d0
info%nan = info%nan/info%nan
info%t_nodal = 1d30

do m = 1,2
	do l = 1,30
		theta = m * (dpd(2)-dpd(1)) + sum(indx(l,m,:)*dpd(:))
		theta = 2d0*theta*rad
		info%ct(l,m) = 2d0*dcos(theta)
		info%st(l,m) = 2d0*dsin(theta)
	enddo
enddo

contains

function csr_read_tide(pathname,uv)
use netcdf
logical :: csr_read_tide
character(*), intent(in) :: pathname
integer(fourbyteint), intent(out) :: uv(:,:,:)
integer(fourbyteint) :: ncid,kx,ky,kw
real(eightbytereal) :: dummy(2)
character(80) :: title

csr_read_tide = .false.

! Open and read netCDF grid of amplitude and phase

if (nf90_open(pathname,nf90_nowrite,ncid) /= nf90_noerr) return
title = ' '
if (nf90_get_att(ncid,nf90_global,"title",title) /= nf90_noerr) return
write (0,'(a)') trim(title)
title = ' '
if (nf90_get_att(ncid,nf90_global,"history",title) /= nf90_noerr) return
write (0,'(a)') trim(title)
if (nf90_inquire_dimension(ncid,1,len=kx) /= nf90_noerr) return
if (nf90_inquire_dimension(ncid,2,len=ky) /= nf90_noerr) return
if (kx /= nx .or. ky /= ny) stop "csrtideinit: Unexpected grid size"
if (nf90_get_att(ncid,1,"actual_range",dummy) /= nf90_noerr) return
if (dummy(1) /= xmin .or. dummy(2) /= xmax) stop "csrtideinit: Unexpected longitude range"
if (nf90_get_att(ncid,2,"actual_range",dummy) /= nf90_noerr) return
if (dummy(1) /= ymin .or. dummy(2) /= ymax) stop "csrtideinit: Unexpected latitude range"
do kw=1,12
   if (nf90_get_var(ncid,kw+2,uv(:,:,kw)) /= nf90_noerr) return
enddo
if (nf90_close(ncid) /= nf90_noerr) return

csr_read_tide = .true.
end function csr_read_tide

end subroutine csrtideinit

!&csrtidefree -- Free up space allocated by csrtideinit
!+
subroutine csrtidefree (info)
type(csrtideinfo), intent(inout) :: info
!
! This routine frees up memory allocated by a previous call to csrtideinit
!
! Argument:
!  info     : Structure initialised by csrtideinit
!-
if (allocated(info%p_ocean)) deallocate (info%p_ocean)
if (allocated(info%p_load)) deallocate (info%p_load)
end subroutine csrtidefree

end module tides_csr
