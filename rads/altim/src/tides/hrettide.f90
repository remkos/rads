module tides_hret
use typesizes
integer(fourbyteint), parameter, private :: nw = 6
real(eightbytereal), parameter, private :: pi = 4d0*atan(1d0), rad = pi/180d0
character(3), parameter, private :: wave(nw) = (/ 'M2 ', 'K1 ', 'S2 ', 'O1 ', 'MA2', 'MB2' /)
type :: hrettideinfo
	integer(fourbyteint) :: nx, ny, undef(nw)
	integer(twobyteint), allocatable :: p_re(:,:,:), p_im(:,:,:)
	real(eightbytereal) :: xmin, xmax, dx, ymin, ymax, dy, fact(nw), arg(nw), f(nw), u(nw), t_nodal, nan, wmin
endtype

contains

!*hrettide -- Compute internal tides according to HRET model
!+
subroutine hrettide (info, utc, lat, lon, tide, tide_comp)
use tides_aux
type(hrettideinfo), intent(inout) :: info
real(eightbytereal), intent(in) :: utc, lat, lon
real(eightbytereal), intent(out) :: tide
real(eightbytereal), intent(out), optional :: tide_comp(nw)
!
! Compute internal tidal height for given time and location from grids
! of harmonic coefficients of the HRET model. The present version
! uses 6 constituents (M2, K1, S2, O1, MA2, MB2).
! This routine is heavily based on the routine PERTH3 by Richard Ray.
!
! The input files can be found in $ALTIM/data/<name>, where $ALTIM
! is an environment variable and <name> is the argument in the
! call to hrettideinit. The grids are in NetCDF format.
!
! To initialize the computation, the function hrettideinit should be
! called first. It allocates the appropriate amount of memory and
! loads the grids into memory. To release the memory for further
! use, call hrettidefree.
!
! Longitude and latitude are to be specified in degrees; time in UTC
! seconds since 1 Jan 1985. All predicted tides are output in meters.
!
! The computation of nodal periods and mean longitudes is optimized
! for the period 1980-2020.
!
! Input arguments:
!  info     : Structure initialised by hrettideinit
!  utc      : UTC time in seconds since 1 Jan 1985
!  lat      : Latitude (degrees)
!  lon      : Longitude (degrees)
!
! Output arguments:
!  tide     : Predicted internal tide using all six constiuents (m)
!  tide_comp: (Optional) Predicted internal tide, separate constituents (m)
!-
! $Log: hrettide.f90,v $
! Revision 1.4  2019/03/28 08:54:50  rads
! - Undefined value (_FillValue) was not initialised or read leading to possible NaN output
!
! Revision 1.2  2019/03/08 09:57:29  rads
! - Introduced internal tide model HRET
!
! (c) Remko Scharroo - EUMETSAT
!-----------------------------------------------------------------------
real(eightbytereal) :: slon
real(eightbytereal), parameter :: mjd85 = 46066d0
integer(fourbyteint) :: istat
complex(eightbytereal) :: val(nw), w(nw)

! When called with invalid time, reset time reference and return

if (.not.(utc < 1d20)) then
	info%t_nodal = 1d30
	tide_comp = info%nan
	tide = info%nan
	return
endif

! If latitude out of range, bail out

if (isnan(lat) .or. isnan(lon) .or. lat < info%ymin .or. lat > info%ymax) then
	tide_comp = info%nan
	tide = info%nan
	return
endif

! Limit longitude to interval of grids

slon = lon
if (slon > info%xmax) then
	slon = slon - 360d0
else if (slon < info%xmin) then
	slon = slon + 360d0
endif

! Compute basic astronomical mean longitudes

call hret_astron (utc)

! Compute all internal tide components

call hret_interp (slon, lat, val, istat)

if (istat == 0) then
	tide_comp = info%nan
	tide = info%nan
else
	w = expj(-(info%arg+info%u)*rad)
	tide_comp = info%f * dble(w*val)
	tide = sum(tide_comp)
endif

contains

!-----------------------------------------------------------------------
! Compute the basic astronomical mean longitudes  s, h, p, N.

subroutine hret_astron (time)
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

info%arg(1) = t2 + 2d0*h - 2d0*s    ! M2
info%arg(2) = t1 + h + 90d0         ! K1
info%arg(3) = t2                    ! S2
info%arg(4) = t1 + h - 2d0*s - 90d0 ! O1
info%arg(5) = info%arg(1) - h       ! MA2
info%arg(6) = info%arg(1) + h       ! MB2

! Determine nodal corrections f and u (only when more than 1 day is passed since last time)

if (abs(time-info%t_nodal) > 86400d0) then
	info%t_nodal = time
	sinn = sin(w*rad)
	cosn = cos(w*rad)
	sin2n = sin(2*w*rad)
	cos2n = cos(2*w*rad)

	info%f( 1) = 1.000d0 - 0.037d0*cosn                 ! M2
	info%f( 2) = 1.006d0 + 0.115d0*cosn - 0.009d0*cos2n ! K1
	info%f( 4) = 1.009d0 + 0.187d0*cosn - 0.015d0*cos2n ! O1

	info%u( 1) = -2.1d0*sinn               ! M2
	info%u( 2) = -8.9d0*sinn + 0.7d0*sin2n ! K1
	info%u( 4) = 10.8d0*sinn - 1.3d0*sin2n ! O1
endif
end subroutine hret_astron

!-----------------------------------------------------------------------
! Interpolates a value from a grid of data at the desired location.
! Interpolation is bilinear.

subroutine hret_interp (x, y, val, istat)
real(eightbytereal), intent(in) :: x, y
complex(eightbytereal), intent(out) :: val(:)
integer(fourbyteint), intent(out) :: istat
integer(fourbyteint) :: kw, i, i0, ii, j, j0, jj
real(eightbytereal) :: wtot, xij, yij, weight(2,2)

istat = 0

! Compute indices for desired position

xij = (x - info%xmin) / info%dx
i0 = floor(xij)
xij = xij - i0

yij = (y - info%ymin) / info%dy
j0 = min(int(yij),info%ny-2) ! Use int() so we do not drop below zero
yij = yij - j0

! Set corner weights

weight(1,1) = (1 - xij) * (1 - yij)
weight(1,2) = (1 - xij) *      yij
weight(2,1) =      xij  * (1-  yij)
weight(2,2) =      xij  *      yij

! Sum up the weighted values for the corners

wtot = 0d0
val = (0d0,0d0)

do i = 1,2
	ii = modulo(i0 + i - 1, info%nx) + 1
	jloop: do j = 1,2
		jj = j0 + j
! Check if the constituents are available:
		if (any(info%p_re(ii,jj,:) == info%undef)) cycle jloop

		istat = istat+1
		wtot = wtot + weight(i,j)

		forall (kw = 1:nw) val(kw) = val(kw) + complex(info%p_re(ii,jj,kw),info%p_im(ii,jj,kw)) * weight(i,j)
	enddo jloop
enddo

if (istat == 0 .or. wtot <= info%wmin) then
	istat = 0
	return
endif
val = val * info%fact / wtot

end subroutine hret_interp

elemental function expj (x)
real(eightbytereal), intent(in) :: x
complex(eightbytereal) :: expj
expj = dcmplx(cos(x),sin(x))
end function expj

end subroutine hrettide

!&hrettideinit -- Initialize GOT tide model
!+
subroutine hrettideinit (pathname, info)
use netcdf
character(*), intent(in) :: pathname
type(hrettideinfo), intent(inout) :: info
!
! Allocate memory for HRET tide modelling and read grids into memory.
!
! Input argument:
!  name     : Pathname of the NetCDF file containing the HRET tide model
!
! Output argument:
!  info     : Structure initialised by hrettideinit
!-
integer(fourbyteint) :: istat, ncid, varid, kw

! Produce log info

write (0,600) trim(pathname)
600 format ('(Loading HRET tide: ',a,')')

! Open and verify netCDF file

if (nf90_open(pathname,nf90_nowrite,ncid) /= nf90_noerr) stop 'hrettideinit: Unable to open file'
if (nf90_inquire_dimension(ncid,1,len=info%nx) + &
    nf90_get_att(ncid,1,'valid_min',info%xmin) + &
    nf90_get_att(ncid,1,'valid_max',info%xmax) /= nf90_noerr) stop 'hrettideinit: Error in longitude dimensions'
if (nf90_inquire_dimension(ncid,2,len=info%ny) + &
    nf90_get_att(ncid,2,'valid_min',info%ymin) + &
    nf90_get_att(ncid,2,'valid_max',info%ymax) /= nf90_noerr) stop 'hrettideinit: Error in latitude dimensions'

! Allocate memory

if (allocated(info%p_re) .or. allocated(info%p_im)) stop 'hrettideinit: Use hrettidefree first'
allocate (info%p_re(info%nx,info%ny,nw), info%p_im(info%nx,info%ny,nw), stat=istat)
if (istat /= 0) stop 'hrettideinit: Not able to allocate memory'

! Initialisations

info%dx = (info%xmax-info%xmin) / (info%nx-1)
info%dy = (info%ymax-info%ymin) / (info%ny-1)
info%fact = 1d0
info%wmin = 0.5d0
info%undef = -32768
info%f = 1d0
info%u = 0d0
info%t_nodal = 1d30
info%nan = 0d0
info%nan = info%nan / info%nan

! Read individual grids

do kw = 1,nw
	if (nf90_inq_varid(ncid,trim(wave(kw))//'re',varid) /= nf90_noerr) stop 'hrettideinit: Missing wave: '//trim(wave(kw))
	if (nf90_get_var(ncid,varid,info%p_re(:,:,kw)) /= nf90_noerr) stop 'hrettideinit: Missing wave: '//trim(wave(kw))
	if (nf90_inq_varid(ncid,trim(wave(kw))//'im',varid) /= nf90_noerr) stop 'hrettideinit: Missing wave: '//trim(wave(kw))
	if (nf90_get_var(ncid,varid,info%p_im(:,:,kw)) /= nf90_noerr) stop 'hrettideinit: Missing wave: '//trim(wave(kw))
	if (nf90_get_att(ncid,varid,'scale_factor',info%fact(kw)) /= nf90_noerr) info%fact(kw) = 1d0
	if (nf90_get_att(ncid,varid,'_FillValue',info%undef(kw)) /= nf90_noerr) info%undef(kw) = -32768
enddo
if (nf90_close(ncid) /= nf90_noerr) return

end subroutine hrettideinit

!&hrettidefree -- Free up space allocated by hrettideinit
!+
subroutine hrettidefree (info)
type(hrettideinfo), intent(inout) :: info
!
! This routine frees up memory allocated by a previous call to hrettideinit
!
! Argument:
!  info : Struct initialised by hrettideinit
!-
if (allocated(info%p_re)) deallocate (info%p_re)
if (allocated(info%p_im)) deallocate (info%p_im)
end subroutine hrettidefree

end module tides_hret
