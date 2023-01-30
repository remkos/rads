!*gimsubs -- Compute TEC from Global Ionosphere Maps
!+
! These FORTRAN90 routines can be used to calculate the vertical total
! electron content based on the Global Ionosphere Maps (GIM) stored
! as annual netCDF files (converted from their IONEX original).
!
! The routines and a brief description:
! gimtec  : Compute TEC at a given location and time
! giminit : Initialize by specifying the location of the data files
!
! To provide the proper FORTRAN90 interface, users need to add the line
! use gimsubs
! to their code.
!
!-----------------------------------------------------------------------
! FUNCTION gimtec (time, lat, lon, info)
! REAL(eightbytereal), INTENT(in) :: time, lat, lon
! TYPE(giminfo), INTENT(inout) :: info
! REAL(eightbytereal) :: gimtec

! Compute the Total Electron Content (TEC) based on the GPS Ionosphere
! Model maps of TEC. This routine loads the appropriate maps from annual
! netCDF files whose location is contained in the <info> structure initialised
! by the giminit function. Once loaded in memory, this function
! performs a spatio-temporal interpolation in the two TEC maps closest
! to the epoch <time>.
!
! Input to the routine are the satellite geographical latitude <lat>,
! longitude <lon> and the epoch of observation in seconds since 1 Jan
! 1985 <time>. Also the structure <info> returned by giminit has to
! be specified.
!
! The TEC is integrated vertically to GPS satellite altitude (20200 km).
! For lower satellites, the TEC needs to be scaled down. Appropriate scaling
! factors are: 0.856 for 800 km (Geosat, GFO, ERS, Envisat) and 0.925 for
! 1350 km (TOPEX, Jason).
!
! The TEC value can be converted to ionospheric path delay by using the
! relationship:
!   PD = C * TEC / f**2
! where PD is path delay in meters, C is a constant of 0.4025 and f is the
! signal frequency in GHz.
!
! To improve efficiency, this routine loads and then keeps in memory the maps
! closest in time to the epoch <time> until one or two new maps need to be
! loaded. These too are stored in the <info> struct.
!
! When either of the two TEC maps is not available for interpolation,
! a NaN value is returned.
!
! During interpolation this routine takes into account the general progression
! of the TEC west ward parallel to the geomagnetic equator.
!
! Input arguments:
!   time   : Time in seconds since 1 Jan 1985
!   lat    : Geographical latitude in degrees
!   lon    : Geographical longitude in degrees
!   info   : Structure initialised by giminit
!
! Output arguments:
!   gimtec : Total Electron Content up to GPS satellite altitude.
!
! Example:
!  pd = 0.925d0 * 0.4025d0 * gimtec (info, time, lat, lon) / 13.6d0**2
!
!-----------------------------------------------------------------------
! FUNCTION gimtec (time, lat, lon)
! CHARACTER(len=*), INTENT(in) :: path
! INTEGER, INTENT(in), OPTIONAL :: verbose
! TYPE(giminfo) :: giminit
!
! This function initialises the use of Global Ionosphere Maps in netCDF format
! for the interpolation of the total electron content contained in those
! files. The netCDF files are yearly files based on the original ascii files
! in the IONEX format produced by JPL, CODE, or other institudes.
!
! The annual files have path names in the form /where/ever/jplg_YYYY.nc, where
! YYYY is the four-digit year number. To use these files, the argument <path>
! to this function must be specified as '/where/ever/jplg_'
!
! The optional argument <verbose> speficies the amount of verbosity of this
! routine (output sent to standard output).
!
! The function returns a structure of type(giminfo) that will be used
! in subsequent calls to gimtec.
!
! Arguments:
!  path    : Path name of the netCDF GIM files, without the YYYY.nc ending
!  verbose : 0 = Errors only (default), 1 = Verbose, 2 = Debugging
!
! Return value:
!  giminit : Structure of type(giminfo) to be used by gimtec
!
! Example:
!  use gimsubs
!  type(giminfo) :: info
!  info = giminit ('/where/ever/jplg_', 1)
!-----------------------------------------------------------------------
! $Log: gimsubs.f90,v $
! Revision 1.1  2010/12/22 22:14:43  rads
! - New version in Fortran 90
! - Renamed from gimtec.f
! - Interface module gimsubs
! - Interpolates in netCDF grids rather than IONEX files
!
! Copyright (c) Remko Scharroo, Altimetrics LLC
!-----------------------------------------------------------------------

module gimsubs
use typesizes

type :: giminfo
	character(len=128) :: path
	real(eightbytereal) :: time(2), time13, dtec
	integer(fourbyteint) :: verbose
	integer(twobyteint) :: tecmap(73,73,2)
endtype

contains

!-----------------------------------------------------------------------
! giminit -- Initialise the processing of Global Ionosphere Maps (netCDF)
!-----------------------------------------------------------------------

function giminit (path, verbose)
use netcdf
character(len=*), intent(in) :: path
integer, intent(in), optional :: verbose
type(giminfo) :: giminit
integer(fourbyteint) :: ncid,varid,hours(2),i
character(len=160) filenm

! Set default values
giminit = giminfo ('', 0d0, 0d0, 0d0, 0, 0)

! Set path and verbosity
giminit%path = path
if (present(verbose)) giminit%verbose = verbose

! Construct full path of the netCDF file for 2002

write (filenm, '(a,a)') trim(path),'2002.nc'

! Check if earlier part of 2002 had odd hours
! If so, set time13 to the transition time of 2-Nov-2002, when IONEX files changed
! from being time-tagged on odd hours to even hours.

if (nf90_open (filenm, nf90_nowrite, ncid) == nf90_noerr) then
	if (nf90_inq_varid (ncid, 'time', varid) == nf90_noerr .and. &
		nf90_get_att (ncid, varid, 'actual_range', hours) == nf90_noerr .and. &
		modulo(hours(1),2) == 1) then
		giminit%time13 = 562809600d0	! 2002-11-02
		if (giminit%verbose >= 2) write (*,550)
	endif
	i = nf90_close (ncid)
endif
550 format ('(giminit: Odd hours before 2002-11-02)')
end function giminit

!-----------------------------------------------------------------------
! gimtec -- Interpolate Global Ionosphere Maps (netCDF)
!-----------------------------------------------------------------------

function gimtec (time, lat, lon, info)
real(eightbytereal), intent(in) :: time, lat, lon
type(giminfo), intent(inout) :: info
real(eightbytereal) :: gimtec
integer(fourbyteint) :: kx, ky, i
real(eightbytereal) :: upy, x, y, dlon, dmap, w(2,2), tec(2)
real(eightbytereal), parameter :: rad=atan(1d0)/45d0, lonp=288.63d0, &
	x0=-180d0, x1=180d0, dx=5d0, y0=-90d0, y1=90d0, dy=2.5d0

! Check if the maps corresponding to the time are loaded.
! - info%time = UTC seconds of the first and second map loaded

if (time < info%time(1) .or. time > info%time(2)) then
	call gimread (1)
	call gimread (2)
endif

! If no grid available, return NaN

if (info%tecmap(1,1,1) == -9999 .or. info%tecmap(1,1,2) == -9999) then
	gimtec = 0d0
	gimtec = gimtec / gimtec
	return
endif

! Factor for latitude adjustment

upy = 0.2d0 * sin((lon-lonp)*rad) * cos(lat*rad)

! Interpolate in first and second grid

do i = 1,2

	! Compute grid indices within the map i (after rotation)

	dlon = (time - info%time(i)) / 240d0
	x = modulo(lon + dlon - x0, 360d0) / dx + 1d0
	kx = floor(x)
	x = x - kx
	y = (lat + dlon * upy - y0) / dy + 1d0
	ky = max(1,min(floor(y),72))
	y = y - ky

	! Determine weights

	w(1,:) = 1d0-x
	w(2,:) = x
	w(:,1) = w(:,1) * (1d0-y)
	w(:,2) = w(:,2) * y

	! Interpolate the TEC in map i

	tec(i) = sum(w * info%tecmap(kx:kx+1,ky:ky+1,i))
enddo

! Interpolate in time between the two values

dmap = (time - info%time(1)) / (info%time(2) - info%time(1))
gimtec = (tec(1) * (1d0-dmap) + tec(2) * dmap) * info%dtec

contains

subroutine gimread (map)
use netcdf
integer(fourbyteint), intent(in) :: map
integer(fourbyteint) :: h,d,y,start(3)=1,ncid,varid,hours(2),odd
real(eightbytereal) :: t
character(len=160) :: filenm
integer(fourbyteint), parameter :: h2000 = 473299200/3600 ! 2000-01-01 in hours since 1985-01-01

550 format ('gimtec: ',a,' ',a)
551 format ('(gimtec: Loading h =',i8,' from ',a,')')
552 format ('(gimtec: Grid for h =',i8,' is already loaded)')
553 format ('(gimtec: Copying h =',i8,')')

! Round the time to the right epoch of the grids

t = time + (map - 1) * 7200d0
if (t < info%time13) then
	odd = 1
else
	odd = 0
endif
h = floor((t/3600d0 - odd) / 2d0) * 2d0 + odd

if (nint(info%time(map)/3600d0) == h) then
	! Skip when grid already loaded
	if (info%verbose >= 2) write (*,552) h-h2000
	return
else if (nint(info%time(3-map)/3600d0) == h) then
	! Copy one grid from the other that is already stored
	if (info%verbose >= 2) write (*,553) h-h2000
	info%time(map) = info%time(3-map)
	info%tecmap(:,:,map) = info%tecmap(:,:,3-map)
	return
endif
info%time(map) = h * 3600d0
info%tecmap(1,1,map) = -9999

! Determine year for current time

d = h / 24
y = floor((d + 0.8d0)/365.25d0) + 1985

! Construct full path of the netCDF file

write (filenm, '(a,i4.4,a)') trim(info%path),y,'.nc'
if (info%verbose >= 1) write (*,551) h-h2000,trim(filenm)

! Open netCDF file

if (nf90_open (filenm, nf90_nowrite, ncid) /= nf90_noerr) then
	write (*,550) 'Error opening',trim(filenm)
	return
endif

! Read start and end time

if (nf90_inq_varid (ncid, 'time', varid) + nf90_get_att (ncid, varid, 'actual_range', hours) /= nf90_noerr) then
	write (*,550) 'Error getting info from',trim(filenm)
	return
endif

! Compute index for this grid

h = h - h2000
if (h < hours(1) .or. h > hours(2)) then
	write (*,550) 'Error: time outside limits of',trim(filenm)
	return
endif
start(3) = (h - hours(1) + 3) / 2
h = nf90_get_var (ncid, varid, info%time(map:map), start(3:3))
info%time(map) = (info%time(map) + h2000) * 3600d0

! Read grid

if (nf90_inq_varid (ncid, 'tec', varid) + nf90_get_att (ncid, varid, 'scale_factor', info%dtec) + &
	nf90_get_var (ncid, varid, info%tecmap(:,:,map:map), start) /= nf90_noerr) then
	write (*,550) 'Error reading TEC map from',trim(filenm)
	info%tecmap(1,1,map) = -9999
	return
endif

h = nf90_close (ncid)
end subroutine gimread

end function gimtec

end module gimsubs
