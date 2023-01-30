program orbite2odr

! Convert CNES ORBITE file to Orbital Data Records (ODR).
!
! syntax: orbite2odr [ options ] odrfile < orbite_file
!-
! $Log: orbite2odr.f90,v $
! Revision 1.9  2020/11/27 16:31:05  rads
! - Added support for Sentinel-6
!
! Revision 1.8  2016/03/10 13:59:07  rads
! - Added Sentinel-3
!
! Revision 1.7  2016/01/20 13:23:25  rads
! - Prepare for Jason-3
!
! Revision 1.6  2014/04/22 11:54:22  rads
! - Added HY-2A (Thanks, John L)
!
! Revision 1.5  2012/05/23 09:58:53  rads
! - Added Jason-1
!
! Revision 1.4  2011/11/04 17:45:02  rads
! - Added another CryoSat format
!
! Revision 1.3  2011/10/25 19:25:28  rads
! - Allow reading of CryoSat orbits
!
! Revision 1.2  2010/11/15 18:10:05  rads
! - Expanded to 4-digit ODR file numbering (arc_????.odr)
!
! Revision 1.1  2010/11/12 01:17:14  rads
! - Based on spx2odr.f90
!
! Created by Remko Scharroo - Altimetrics LLC
!***********************************************************************
use typesizes
character(80) :: odrname=' ',arg,altim
character(4) :: spec='xODR'
character(8) :: satel=''
character(36) :: month_string='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'
integer(fourbyteint) :: lstr,i,krec=0,itime0,itime1,iargc, &
		irem=0,irep=0,iarc=0,iskip=0,mjd,ios
real(eightbytereal) :: ecf(3),utc,tai,skip=-1d0,dlon,dlat,dhgt,r
real(eightbytereal), parameter :: pi=4d0*atan(1d0), rad=pi/180d0, mjd85=46066d0
logical :: ltlend

! Variables for ODR format

type :: odr
	integer(fourbyteint) :: itime,lat,lon,hgt
endtype
type(odr) :: o

! Variables for ORBITE format

character(132) :: line
integer(fourbyteint) :: yyyy, month, day, hour, minute
real(eightbytereal) :: interval = 60d0, seconds
character(3) :: mon
logical :: eef = .false.

! Variables for TAI

integer(fourbyteint), parameter :: nutc=100
real(eightbytereal) :: tai_utc(nutc)
integer(fourbyteint) :: mjdutc(nutc)
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
		satel = arg(5:)
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
1301 format ('orbite2odr: Convert CNES ORBITE orbit file to ODR file'// &
'syntax: orbite2odr [options] < orbite_file odrfile'// &
'where [options] are:'/ &
'  sat=SATEL   : specify satellite'/ &
'  rem=REMID   : specify orbit ID'/ &
'  skip=DAY    : specify days to skip at start of orbit')
	stop
endif

! Open ODR file

open (20,file=odrname,status='replace',form='unformatted',recl=16,access='direct')
i=index(odrname,'ODR.')
if (i.gt.0) read (odrname(i+4:i+6),*) iarc
i=index(odrname,'.odr')
if (i.gt.0) read (odrname(i-4:i-1),*) iarc
call setearth(6378136.3d0,298.257d0)

! Read ORBITE header records

read (*,'(a)') line
if (line(:2) == '20') then
	! This is a headerless CRYOSAT file
	eef = .true.
	satel = 'CRYOSAT2'
else
	! Parse the product headers
	do
		read (*,'(a)') line
		if (line(:14) == 'SPH_DESCRIPTOR') then
			if (index(line,'ENV1') > 0) then
				satel='ENVISAT1'
				irep=35000
			else if (index(line,'jason1') > 0) then
				satel='JASON-1'
				irep=9920
			else if (index(line,'jason2') > 0) then
				satel='JASON-2'
				irep=9920
			else if (index(line,'jason3') > 0) then
				satel='JASON-3'
				irep=9920
			else if (index(line,'cryosat') > 0) then
				satel='CRYOSAT2'
			else if (index(line,'hy2a') > 0) then
				satel='HY-2A'
				irep=14000
			else if (index(line,'sent3a') > 0) then
				satel='SNTNL-3A'
				irep=27000
			else if (index(line,'sent3b') > 0) then
				satel='SNTNL-3B'
				irep=27000
			else if (index(line,'sent6a') > 0) then
				satel='SNTNL-6A'
				irep=9920
			else if (index(line,'sent6b') > 0) then
				satel='SNTNL-6B'
				irep=9920
			endif
		else if (line == 'DSR_SIZE=+0000000129<bytes>') then
			read (*,*)
			exit
		else if (line == 'DSR_SIZE=+0000000097<bytes>') then
			read (*,*)
			read (*,*)
			exit
		endif
	enddo
	read (*,'(a)') line
endif
if (satel == '') stop 'Unknown satellite'

! Read ORBITE epoch and position/velocity records
! Note that we have one line still in memory

do
	if (eef) then
		read (line,'(i4,1x,i2,1x,i2,i3,1x,i2,1x,f9.6,7x,3f13.3)') yyyy,month,day,hour,minute,seconds,ecf
		call ymd2mjd (yyyy,month,day,mjd)
		utc = (mjd-mjd85)*86400d0+hour*3600d0+minute*60d0+seconds
	else if (satel == 'ENVISAT1') then
		read (line,'(i2,1x,a3,1x,i4,i3,1x,i2,1x,f9.6,16x,3f13.3)') day,mon,yyyy,hour,minute,seconds,ecf
		month = index(month_string,mon) / 3 + 1
		call ymd2mjd (yyyy,month,day,mjd)
		utc = (mjd-mjd85)*86400d0+hour*3600d0+minute*60d0+seconds
	else
		read (line,'(i5,f12.6,3f13.4)') day,seconds,ecf
		tai = (day - 12784) * 86400d0 + seconds
		utc = tai - tai_utc(iutc)
		! Convert TAI to UTC
		do while (utc/86400d0 >= mjdutc(iutc+1) .and. modulo(utc,86400d0) > 0.5d0 .and. iutc < nutc)
			iutc = iutc + 1
		enddo
		utc = tai - tai_utc(iutc)
	endif

	o%itime = nint(utc)

	call xyzgeo(ecf,r,dlat,dlon,dhgt)
	dlat=dlat/rad
	dlon=dlon/rad
	if (dlon > 180d0) dlon=dlon-360d0
	krec=krec+1

	o%lat=nint(dlat*1d7)
	o%lon=nint(dlon*1d7)
	o%hgt=nint(dhgt*1d3)
	write (20,rec=krec+2) o

	read (*,'(a)',iostat=ios) line
	if (ios /= 0 .or. index(line,'EOF') > 0) exit
enddo

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

end program orbite2odr
