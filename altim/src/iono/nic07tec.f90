!*NIC07TEC -- Evaluate NIC at given location and time
!+
function nic07tec (time, lat, lon, nicfile, fluxfile)
use typesizes
use netcdf
real(eightbytereal), intent(in) :: time,lat,lon
character(*), intent(in) :: nicfile,fluxfile
real(eightbytereal) :: nic07tec
!
! This function determines the total electron content at a given time (time)
! and location (lat, lon) based on the NOAA Ionosphere Climatology (NIC).
!
! Input parameters are:
!  time    : time in UTC seconds since 1.0 Jan 1985
!  lat     : latitude in degrees
!  lon     : longitude in degrees (range from -540 to 540 is allowed)
!  nicfile : path name of the NIC climatology file
!  fluxfile: path name of the radio flux file
!
! Returned value:
!  nic07tec: total electron content (TEC units = 10**16 electrons/m**2)
!
! The two input files are loaded on the first call of the function. On any
! subsequent call the last two arguments are ignored.
!
! The routine will halt the processing when the requested time is out of
! reach of the flux file. Hence, the flux file should be regularly updated.
! To allow some extrapolation, 7 days of flux values are automatically added
! to the end of the flux file.
!
! The TEC is integrated to the height mentioned in the climatology files.
! Generally, two climatology files will be available (e.g. nic07_heo.nc and
! nic07_leo.nc) for high and low earth orbiter, respectively.
!
! The climatology files and flux file can be found at
! ftp://falcon.grdl.noaa.gov/pub/remko/nic07
!
! The TEC value can be converted to ionospheric path delay by using the
! relationship:
!   PD = C * TEC / f**2
! where PD is path delay in metres, C is a constant of 40250d13 and f is the
! signal frequency in Hz.
!-
!  $Log: nic07tec.f90,v $
!  Revision 1.1  2008/01/22 17:51:40  rads
!  - Renamed nictec to nic07tec
!
!  Revision 1.1  2007/04/06 13:56:31  rads
!  - Initial version
!
!-----------------------------------------------------------------------
integer, parameter :: md=12000
real(eightbytereal), save :: flux(md),tec_scale,flux_range(2)
integer(twobyteint), save :: nic(73,73,12,12,2)
integer, save :: nd=-1
integer :: i,ncid,varid,kd,h0,h1,m0,m1,kx,ky
real(eightbytereal) :: xflux,xtime,d,h,m,x,y,tec_cube(2,2,2)
real(eightbytereal), save :: day0=0d0
logical :: isopen

if (nd < 0) then	! Initialisation

! Read NIC climatology

   call nfs(nf90_open(nicfile,nf90_nowrite,ncid))
   call nfs(nf90_inq_varid(ncid,"flux",varid))
   call nfs(nf90_get_var(ncid,varid,flux_range))
   call nfs(nf90_inq_varid(ncid,"tec",varid))
   call nfs(nf90_get_att(ncid,varid,"scale_factor",tec_scale))
   call nfs(nf90_get_var(ncid,varid,nic))
   call nfs(nf90_close(ncid))

! Read flux data

   do ncid=99,7,-1
      inquire (unit=ncid,opened=isopen)
      if (.not.isopen) exit
   enddo
   open (ncid,file=fluxfile)
   nd = 0
   do
      read (ncid,*,iostat=i) xtime,xflux
      if (i /= 0) exit
      nd = nd + 1
      if (nd > md) call fin ("Too many values in flux file")
      if (nd == 1) day0 = xtime
      flux(nd) = xflux
   enddo
   close (ncid)

! Copy last record 7 times to allow some extrapolation

   flux(nd+1:nd+7) = flux(nd)

endif

! Convert time (UTC85) into days from 2000.0

xtime = (time - 473299200d0)/86400d0

! Determine day fraction

d = xtime - day0 + 1d0
kd = floor(d)
d = d - kd
if (kd < 1 .or. kd > nd + 6) call fin ("Time outside flux file")

! Interpolate the flux and determine flux fraction

xflux = flux(kd) * (1d0 - d) + flux(kd+1) * d
xflux = (xflux - flux_range(1)) / (flux_range(2) - flux_range(1))

! Determine month and month fraction. m0 and m1 are neighbouring months.

m = modulo (xtime / 30.4375d0 + 0.5d0, 12d0)
m0 = floor(m)
m = m - m0
m1 = m0 + 1
if (m0 == 0) m0 = 12

! Determine hour and hour fraction. h0 and h1 are neighbouring hours.

h = modulo (xtime, 1d0) * 12d0 + 1
h0 = floor(h)
h = h - h0
h1 = h0 + 1
if (h1 == 13) h1 = 1

! Determine Y coordinate and Y fraction

y = (lat + 90d0)/2.5d0 + 1
ky = floor(y)
y = y - ky

! Determine X coordinate, shifted East for earlier grid, and X fraction

x = lon + h*30d0
if (x < -180d0) x = x + 360d0
if (x >= 180d0) x = x - 360d0
x = (x + 180d0)/5d0 + 1
kx = floor(x)
x = x - kx

! Store earlier grid data while interpolating in month and weighting
! by hour, so tec_cube will be 2x2x2 matrix of which corners are the
! neighbouring longitude nodes, latitude nodes and low-high flux level.

tec_cube(:,:,:) = (1-h) * (nic(kx:kx+1,ky:ky+1,h0,m0,:)*(1-m) &
		+ nic(kx:kx+1,ky:ky+1,h0,m1,:)*m)

! Determine X coordinate, shifted West for later grid (fraction stays the same).

kx = kx - 6
if (kx < 1) kx = kx + 72

! Add later grid data while interpolating in month and weighting by hour

tec_cube(:,:,:) = tec_cube(:,:,:) + h * (nic(kx:kx+1,ky:ky+1,h1,m0,:)*(1-m) &
		+ nic(kx:kx+1,ky:ky+1,h1,m1,:)*m)

! Interpolate in the remaining dimensions (longitude, latitude and flux)

tec_cube(1,:,:) = (1-x)*tec_cube(1,:,:) + x*tec_cube(2,:,:)
tec_cube(1,1,:) = (1-y)*tec_cube(1,1,:) + y*tec_cube(1,2,:)
nic07tec = ((1-xflux)*tec_cube(1,1,1) + xflux*tec_cube(1,1,2)) * tec_scale

contains

subroutine nfs(ios)
integer, intent(in) :: ios
if (ios==nf90_noerr) return
call fin(nf90_strerror(ios))
end subroutine nfs

subroutine fin(string)
character(*), intent(in) :: string
character(80) :: prognm
call getarg(0,prognm)
write (0,'(a,": nic07tec: ",a)') trim(prognm),trim(string)
stop
end subroutine fin

end function nic07tec
