!*nicsubs -- Evaluate NIC at given location and time
!+
! These FORTRAN90 routines can be used to calculate the vertical total
! electron content based on the NOAA Ionosphere Climatology (NIC).
! The versions NIC08 and NIC09 (and possibly later versions) are
! supported.
!
! The routines and a brief description:
! nictec  : Evaluate NIC at a given location and time
! nicinit : Initialize by reading the climatology and forcing files
! nicfree : Free up allocated memory
!
! To provide the proper FORTRAN90 interface, users need to add the line
! use nicsubs
! to their code.
!
! Reference:
! Scharroo, R., and W. H. F. Smith, A global positioning system-based
! climatology for the total electron content in the ionosphere,
! J. Geophys. Res., 115, A10318, 2010, doi:10.1029/2009JA014719.
!
!-----------------------------------------------------------------------
! FUNCTION nictec (time, lat, lon)
! REAL(eightbytereal), INTENT(in) :: time, lat, lon
! REAL(eightbytereal) :: nictec
!
! This function determines the vertical total electron content at a given
! time (time) and location (lat, lon) based on the NOAA Ionosphere Climatology
! (NIC).
!
! The nictec routine needs to be initilized using a call to nicinit. This
! routine will load the NIC climatology file (e.g., nic09_clim.nc) and a
! regularly updated file (e.g., nic09_gtec.nc) containing the forcing
! parameter for this climatology: 2-hourly values of the global mean TEC.
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
! Input arguments are:
!  time    : time (UTC seconds since 1.0 Jan 1985)
!  lat     : latitude (degrees) 
!  lon     : longitude (degrees) (range from -540 to 540 is allowed)
!
! Returned value:
!  nictec  : total electron content (TEC units = 10^16 electrons/m^2)
!
! Example:
!  pd = 0.925d0 * 0.4025d0 * nictec (time, lat, lon) / 13.6d0**2
!
!-----------------------------------------------------------------------
! SUBROUTINE nicinit (nicfile, gtecfile)
! CHARACTER(len=*), INTENT(in) :: nicfile, gtecfile
!
! This routine needs to be called prior to nictec. It reads the climatology
! file and the forcing parameter file into memory and performs some other
! initializations.
!
! The climatology files and forcing parameter file can be found at
! ftp://falcon.grdl.noaa.gov/pub/remko/nic09 (latest version)
! ftp://falcon.grdl.noaa.gov/pub/remko/nic08 (older version)
!
! The nictec routine will return NaN (not-a-number) when the requested time
! is out of reach of the forcing parameter file (gtecfile). Hence that files
! should be updated reguraly. You can find an updated version of this file
! on the abovementioned FTP site daily.
!
! When nictec is no longer needed in your program, the memory allocated
! by nicinit can be freed using the nicfree routine.
!
! Input arguments are:
!  nicfile : path name of the NIC climatology file
!  gtecfile: path name of the GTEC forcing parameter time series
!            or use 'gtec=<value>' to set a constant value.
!
! Example:
!  use nicsubs
!  call nicinit ('nic09_clim.nc', 'nic09_gtec.nc')
!
!-----------------------------------------------------------------------
! SUBROUTINE nicfree ()
!
! Free up the memory previously allocated by nicinit.
! The nicfree routine has no arguments.
!
!-----------------------------------------------------------------------
! MODULE nicsubs
!
! This module must be added to the calling routine by entering the line
! "use nicsubs" in your code. It gives also access to NIC's internal
! variables, which are all contained in the structure
! nic of type nicinfo, with the elements:
!
! tec_scale   : Scale factor to TECU for climatology
! gtec_scale  : Scale factor to TECU for GTEC
! level       : Lower and upper level of GTEC and difference
! trange      : Time range in the GTEC table (hours since 2000.0) and time step
! gtec        : GTEC value at current time (after nictec call)
! dtec        : Derivative of TEC to GTEC (after nictec call)
! nt          : Number of entries in the GTEC table (0 when uninitialized)
! year        : First and last year of JPL GIM data used in NIC
! gtec_tab    : Table of global mean TEC forcing
! tec         : Array of the NIC climatology data
! 
!-----------------------------------------------------------------------
! $Log: nicsubs.f90,v $
! Revision 1.6  2010/12/22 15:01:49  rads
! - New interface (module nicsubs)
! - New interpolation scheme, fits better with NIC development
!    (produces only minor changes with previous version, level 0.1 TECU)
!
! Revision 1.4  2009/04/20 16:35:27  remko
! - Improved description (aimed at NIC09)
! - Return NaN when time out of range
! - Otherwise compatible with previous version
!
! Revision 1.3  2008/01/22 17:48:15  remko
! - Use more unique variable names in module
!
! Revision 1.2  2008/01/15 20:41:26  remko
! - New version for NIC08. Includes nicinit and nicfree.
!
! Revision 1.1  2007/04/20 14:52:57  remko
! - Initial revision
!
! Copyright (c) Remko Scharroo, Altimetrics LLC
!-----------------------------------------------------------------------

module nicsubs
use typesizes
type :: nicinfo
	real(eightbytereal) :: tec_scale,gtec_scale,level(3),trange(3),gtec,dtec
	integer(fourbyteint) :: nt,year(2)
	integer(twobyteint), allocatable :: gtec_tab(:),tec(:,:,:,:,:)
endtype
type(nicinfo) :: nic = nicinfo(0d0,0d0,0d0,0d0,0d0,0d0,0,0,null(),null())

private :: nicfin

contains

!-----------------------------------------------------------------------
! nicinit: Load NIC climatology and GTEC forcing file
!-----------------------------------------------------------------------

subroutine nicinit (nicfile, gtecfile)
use netcdf
character(len=*), intent(in) :: nicfile,gtecfile
integer(fourbyteint) :: ncid,varid
character(len=32) :: text
real(eightbytereal) :: gtec

if (nic%nt > 0) call nicfin('nicinit: NIC already initialized. Use nicfree first.')

! Read NIC climatology

call nfs(nf90_open(nicfile,nf90_nowrite,ncid))
call nfs(nf90_get_att(ncid,nf90_global,'description',text))
if (text(:5)=='NIC07') call nicfin('nicinit: This version does not work with NIC07')
call nfs(nf90_get_att(ncid,nf90_global,'year_range',nic%year))
call nfs(nf90_inq_varid(ncid,'gtec',varid))
call nfs(nf90_get_var(ncid,varid,nic%level(1:2)))
call nfs(nf90_inq_varid(ncid,'tec',varid))
call nfs(nf90_get_att(ncid,varid,'scale_factor',nic%tec_scale))
allocate (nic%tec(73,73,12,12,2))
call nfs(nf90_get_var(ncid,varid,nic%tec))
call nfs(nf90_close(ncid))

! Make nic%level(3) the GTEC range (to save some computation later on)

nic%level(3) = nic%level(2) - nic%level(1)

! If gtecfile starts with 'gtec=', simply set a constant GTEC

if (gtecfile(:5) == 'gtec=') then
    nic%nt=2
    nic%gtec_scale = 1d-2
    nic%trange = (/-1d20,2d20,0d0/)
    allocate (nic%gtec_tab(nic%nt))
    read (gtecfile(6:),*) gtec
    nic%gtec_tab = nint(gtec/nic%gtec_scale)
    return
endif

! Read GTEC data

call nfs(nf90_open(gtecfile,nf90_nowrite,ncid))
call nfs(nf90_inquire_dimension(ncid,1,len=nic%nt))
allocate (nic%gtec_tab(nic%nt))
call nfs(nf90_inq_varid(ncid,'time',varid))
call nfs(nf90_get_att(ncid,varid,'actual_range',nic%trange(1:2)))
call nfs(nf90_inq_varid(ncid,'gtec',varid))
call nfs(nf90_get_att(ncid,varid,'scale_factor',nic%gtec_scale))
call nfs(nf90_get_var(ncid,varid,nic%gtec_tab))
call nfs(nf90_close(ncid))

! Make nic%trange(3) the time step (to save some computation later on)

nic%trange(3) = (nic%trange(2) - nic%trange(1)) / (nic%nt - 1)

contains

subroutine nfs(ios)
! Error exit for NetCDF calls
integer, intent(in) :: ios
if (ios == nf90_noerr) return
call nicfin('nicinit: '//nf90_strerror(ios))
end subroutine nfs

end subroutine nicinit

!-----------------------------------------------------------------------
! nicfin: General exit for all NIC routines
!-----------------------------------------------------------------------

subroutine nicfin(string)
character(len=*), intent(in) :: string
character(len=80) :: prognm
call getarg(0,prognm)
write (0,'(a,": ",a)') trim(prognm),trim(string)
stop
end subroutine nicfin

!-----------------------------------------------------------------------
! nicfree: Free up memory claimed by nicinit
!-----------------------------------------------------------------------

subroutine nicfree ()
if (nic%nt > 0) deallocate (nic%gtec_tab,nic%tec)
nic%nt = 0
end subroutine nicfree

!-----------------------------------------------------------------------
! nictec: Interpolate NIC at given time, latitude and longitude
!-----------------------------------------------------------------------

function nictec (time, lat, lon)
real(eightbytereal) :: nictec
real(eightbytereal), intent(in) :: time, lat, lon
integer(fourbyteint) :: kt
real(eightbytereal) :: t, xtime, upy
real(eightbytereal), parameter :: rad=atan(1d0)/45d0, lonp=288.63d0

if (nic%nt <= 0) call nicfin('nictec: NIC not initialized. Use nicinit first.')

! Convert time (UTC85) into hours since 2000.0

xtime = (time - 473299200d0) / 3600d0

! Determine location of matching GTEC record

t = (xtime - nic%trange(1)) / nic%trange(3) + 1d0
kt = floor(t)

! When the time is out of range, return NaN

if (kt < 1 .or. kt >= nic%nt) then
    nictec = 0d0
    nictec = nictec / nictec
    return
endif

upy = 0.2d0 * sin((lon-lonp)*rad) * cos(lat*rad)

! Interpolate at the times kt and kt+1

nic%gtec = 0d0
nic%dtec = 0d0
nictec = 0d0
call nictec_inter (kt)
call nictec_inter (kt+1)

contains

subroutine nictec_inter (kt)
integer(fourbyteint), intent(in) :: kt
integer(fourbyteint) :: h0,m0,m1,kx,ky,j
real(eightbytereal) :: tec(2), dmap, gtec, f, dlon, w(2,2), x, y, m, xtime
real(eightbytereal), parameter :: x0=-180d0, x1=180d0, dx=5d0, y0=-90d0, y1=90d0, dy=2.5d0

! Interpolate forcing parameter (GTEC) and determine its fraction.

gtec = nic%gtec_tab(kt) * nic%gtec_scale
f = (gtec - nic%level(1)) / nic%level(3)
xtime = nic%trange(1) + (kt-1) * nic%trange(3)

! Determine month and month fraction. m0 and m1 are neighbouring months.

m = modulo (xtime / 730.5d0 + 0.5d0, 12d0)
m0 = floor(m)
m = m - m0
m1 = m0 + 1
if (m0 == 0) m0 = 12

! Determine bihour index

h0 = modulo(kt-1,12)+1

! Compute grid indices within the TEC map, including longitude shift

dlon = (t - kt) * 30d0
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

! Store grid data while interpolating in latitude, longitude and month,
! so that tec_sqr will be a 2x2 matrix of which corners are the
! neighbouring bihours and forcing level.

forall (j = 1:2)
	tec(j) = (1d0-m) * sum(w * nic%tec(kx:kx+1,ky:ky+1,h0,m0,j)) + &
	              m  * sum(w * nic%tec(kx:kx+1,ky:ky+1,h0,m1,j))
endforall
dmap = 1d0 - abs(t-kt)
nic%gtec = nic%gtec + dmap * gtec
nic%dtec = nic%dtec + dmap * (tec(2) - tec(1)) * nic%tec_scale / nic%level(3)
nictec = nictec + dmap * ((1d0-f) * tec(1) + f * tec(2)) * nic%tec_scale
end subroutine nictec_inter

end function nictec

end module nicsubs
