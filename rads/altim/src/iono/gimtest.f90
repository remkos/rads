program gimtest

! Program to test the GIM software.
!
! Output are the TEC values for epochs between 1 Nov 2002 00:00 and 3 Nov 2002 0:00
! at location 31.1 N, 18.6 E.
!-----------------------------------------------------------------------
use typesizes
use gimsubs
use nicsubs

integer(fourbyteint) :: i
real(eightbytereal) :: tec(4),mjd=52579d0,lat=31.1d0,lon=18.6d0,f,freq=13.575d0,utc,gimtec_old
type(giminfo) :: info(2)

f=-0.40250d0/freq**2

info(1) = giminit('/rads/altim/data/gim/jplg_',2)
info(2) = giminit('/rads/altim/data/nic09/nic09_',2)
call nicinit ('/rads/altim/data/nic09/nic09_clim.nc','/rads/altim/data/nic09/nic09_gtec.nc')

do i=0,48*60,30
	utc=(mjd-46066)*86400d0+i*60
	tec(1)=gimtec_old(utc,lat,lon)
	tec(2)=gimtec(utc,lat,lon,info(1))
	tec(3)=nictec(utc,lat,lon)
	tec(4)=gimtec(utc,lat,lon,info(2))
	write (*,'(f12.1,8f9.4)') utc,tec,f*tec
enddo

end program gimtest

function gimtec_old(utc,lat,lon)
use typesizes
real(eightbytereal) :: utc,lat,lon,gimtec_old
real(eightbytereal) :: gimtec
integer(fourbyteint), parameter :: ndays=2
integer(twobyteint) :: tecmap(73,71,ndays*12+1)
gimtec_old = gimtec(utc,lat,lon,ndays,tecmap,.true.,'JPLG')
end function gimtec_old
