program prc2odr
!
! This program converts the GFZ PRC orbit format to ODRs
!
! syntax: prc2odr odrfile < prcfile
!-
! 17-Aug-17 - Remko Scharroo (DUT/DEOS)
! $Log: prc2odr.f90,v $
! Revision 1.2  2012/06/11 19:14:36  rads
! - Open ODR file as "replace" not "new"
!
! Revision 1.1  2010/06/14 19:18:24  rads
! - Revamped version of prc2odr.f
!
!-----------------------------------------------------------------------
use typesizes
character(80) :: odrname=' ',arg,altim
character(132) :: line
character(4) :: spec='xODR'
character(8) :: satel='ERS-1'
integer(fourbyteint) :: lstr,i,krec=0,irep=35000,iarc,irem=9055,isat
integer(fourbyteint) :: itime0,itime1,iskip,ios
real(eightbytereal) :: ecf(6),et,utc,day,sec,dlat,dlon,r,dhgt,interval=60d0,skip=-1d0
real(eightbytereal),parameter :: pi=4d0*atan(1d0), rad=pi/180d0, masrad=rad/3600d3, &
	xpole0=45d0,ypole0=286d0,mjd2000=51544.5d0-46066d0
logical :: ltlend

! Variables for ODR format

type :: odr
	integer(fourbyteint) :: itime,lat,lon,hgt
endtype
type(odr) :: o

! Variables for TAI

integer(fourbyteint), parameter :: nutc=100
real(eightbytereal) :: tai_utc(nutc)
integer(fourbyteint) :: mjdutc(nutc)
real(eightbytereal), parameter :: et_tai = 32.184d0
integer(fourbyteint) :: iutc

! Read TAI - UTC table

call checkenv('ALTIM',altim,lstr)
open (10,status='old',form='formatted',file=altim(1:lstr) // '/data/tables/tai-utc.dat')
iutc=1
do
	read (10,'(16x,i8,12x,f12.7)',iostat=ios) mjdutc(iutc),tai_utc(iutc)
	if (ios /= 0) exit
	mjdutc(iutc)=mjdutc(iutc)-2446066
	iutc=iutc+1
enddo
close (10)
mjdutc(iutc)=999999
tai_utc(iutc)=tai_utc(iutc-1)
iutc=1

! Scan the command line

do i=1,iargc()
	call getarg(i,arg)
	if (arg(:4) == 'sat=') then
		satel=arg(5:)
	else if (arg(:4) == 'rem=') then
		read (arg(5:),*) irem
	else if (arg(:5) == 'skip=') then
		read (arg(6:),*) skip
	else
		odrname=arg
	endif
enddo

! Write usage info on faulty input

if (odrname == ' ') then
	write (*,1301)
1301 format ('prc2odr: Convert PRC orbit file to ODR file'// &
'syntax: prc2odr [options] < prcfile odrfile'// &
'where [options] are:'/ &
'  rem=REMID   : specify orbit ID'/ &
'  skip=DAY    : specify days to skip at start of orbit')
	stop
endif

call setearth(6378136.3d0,298.257d0)

! Open ODR file

open (20,file=odrname,status='replace',form='unformatted',recl=16,access='direct')
i=index(odrname,'ODR.')
read (odrname(i+4:i+6),*) iarc

! Convert PRC file

do
	read (*,550,iostat=ios) line
	if (ios /= 0) exit
	if (line(1:6) /= 'STTERR') cycle
	read (line,20) isat,day,sec,ecf

! Convert Ephemeris time to UTC
! UTC - ET = (TAI - ET) - (TAI - UTC)

	et = (day + mjd2000)*86400d0 + sec
	utc = et - et_tai - tai_utc(iutc)

	do while (utc/86400d0 >= mjdutc(iutc+1) .and. modulo(utc,86400d0) > 0.5d0 .and. iutc < nutc)
		iutc=iutc+1
	enddo
	utc = utc - et_tai - tai_utc(iutc)
	
	if (krec == 1) interval = utc - o%itime
	o%itime = nint(utc)
	
	call rotate(2,-xpole0*masrad,ecf,ecf)
	call rotate(1,-ypole0*masrad,ecf,ecf)
	
	call xyzgeo(ecf,r,dlat,dlon,dhgt)
	dlat=dlat/rad
	dlon=dlon/rad
	if (dlon > 180d0) dlon=dlon-360d0
	krec=krec+1
	
	o%lat=nint(dlat/rad*1d7)
	o%lon=nint(dlon/rad*1d7)
	o%hgt=nint(dhgt*1d3)
	write (20,rec=krec+2) o
enddo

20  format (6x,i7,1x,f6.1,f11.6,3f12.3,3f11.3)
550 format (a)

! Reread the first and last two records.

read (20,rec=3) itime0
read (20,rec=krec+1) itime1
read (20,rec=krec+2) o

! Check if there was a leapsecond in the data period.
! If it is on the last record only, do not apply it there.

if (modulo(o%itime-itime1,10) == 9) then
	o%itime=o%itime+1
	write (20,rec=krec+2) o
endif

! Search for begin of precise arc.
! If skip = 0: start immediately
! If skip < 0: start immediately only when arc does not start at noon or midnight
! Else: start around abs(skip) from beginning or arc

if (skip == 0d0) then
	iskip=itime0
else if (modulo(itime0,43200) /= 0 .and. skip < 0d0) then
	iskip=itime0
else
	skip=itime0+abs(skip)*86400
	call search(dble(itime0),interval,skip,krec)
	iskip=nint(skip)
endif

! If little endian: read the whole file again and swap the integers

if (ltlend()) then
	do i=3,krec+2
		read (20,rec=i) o
		call i4swap(4,o)
		write (20,rec=i) o
	enddo
	call i4swap(1,iskip)
	call i4swap(1,irep)
	call i4swap(1,iarc)
	call i4swap(1,krec)
	call i4swap(1,irem)
endif

! Determine satellite

select case (isat)
case (9105001)
	satel='ERS-1'
case default
	satel='ERS-2'
end select

! Write headers

write (20,rec=1) spec,satel,iskip
write (20,rec=2) irep,iarc,krec,irem
close (20)

contains

subroutine search(time0,step,skip,krec)
real(eightbytereal) time0,step,skip,dlat,x
integer(fourbyteint) :: krec,irec,itime0,itime1,itime2,lat0,lat1,lat2

irec=min(nint((skip-time0)/step)+1,krec-1)
do
	read (20,rec=2+irec-1) itime0,lat0
	read (20,rec=2+irec  ) itime1,lat1
	read (20,rec=2+irec+1) itime2,lat2
	if (lat1 < lat0 .and. lat1 < lat2 .or. irec == krec-1) then
		call inter(dble(lat0),dble(lat1),dble(lat2),x,dlat)
		skip=itime1+x*step
		return
	else if (lat2 < lat0) then
		irec=irec+1
	else
		irec=irec-1
	endif
enddo
end subroutine search

subroutine inter(x0,x1,x2,t,xt)
real(eightbytereal) :: x0,x1,x2,t,xt,b,c
b=(x2-x0)/2
c=x2-b-x1
t=-b/(2*c)
xt=x1+b*t+c*t*t
end subroutine inter

end program prc2odr
