program poe2odr

! Convert NASA POE file to Orbital Data Records (ODR).
!
! syntax: poe2odr < poefile odrfile
!
! A sample POE record is:
!
! 0.2001060100000000D+120.0000000000000000D+000.2495428396768284D+030.2243999999999997D+030.3939999999999995D+030.1520003707683067D+03
! 0.6295801307003660D+07-.6521191703429127D+06-.3367460041527528D+070.2989305324217648D+04-.2940451954818308D+040.6163898617890231D+04
! -.1589433772339607D+070.6126669061810342D+07-.3367460041527528D+070.2156989600723170D+040.3944392839016151D+040.6163898617890231D+04
! 00000000000000000000000.0000000000000000D+000.0000000000000000D+000.0000000000000000D+000.0000000000000000D+00
!
!-
! $Log: poe2odr.f90,v $
! Revision 1.1  2012/06/11 19:13:43  rads
! - Converted poe2odr.f to Fortran 90
!
! Created by Remko Scharroo - Altimetrics LLC
!-----------------------------------------------------------------------
use typesizes
character(80) :: odrname=' ',arg
character(4) :: spec='xODR'
character(8) :: satel='GFO-1'
integer(fourbyteint) :: i,krec=0,itime0,itime1,iargc,irem=0,irep,iarc,iskip,mjd,ios,interval,leapsec,mdate
real(eightbytereal) :: ecf(3),xyz(3),skip=-1d0,dlon,dlat,dhgt,r,begin=0d0,sec85,xpole,ypole,grhran,sec,time
integer(fourbyteint), parameter :: mjd85=46066
real(eightbytereal), parameter :: pi=4d0*atan(1d0), rad=pi/180d0, masrad=rad/3600d3
logical :: ltlend

! Variables for ODR format

type :: odr
	integer(fourbyteint) :: itime,lat,lon,hgt
endtype
type(odr) :: o

! Scan the command line

do i=1,iargc()
	call getarg(i,arg)
	if (arg(:4) == 'sat=') then
		satel=arg(5:)
	else if (arg(:4) == 'rem=') then
		read (arg(5:),*) irem
	else if (arg(:5) == 'skip=') then
		read (arg(6:),*) skip
	else if (arg(:6) == 'begin=') then
		read (arg(7:),*) begin
		begin = sec85(0,begin)
	else
		odrname=arg
	endif
enddo

! Write usage info on faulty input

if (odrname == ' ') then
	write (*,1301)
1301 format ('poe2odr: Convert NASA POE file to ODR file'// &
'syntax: poe2odr [options] < poefile odrfile'// &
'where [options] are:'/ &
'  sat=SATEL   : specify 8-character satellite name (def: GFO-1)'/ &
'  rem=REMID   : specify orbit type (def: 777)')
	stop
endif

! Check satellite

if (satel(:3) == 'GFO' .or. satel(:3) == 'GEO') then
   irep=17050
else if (satel(:3) == 'TOP' .or. satel(:3) == 'JAS') then
   irep=9920
else if (satel(:3) == 'ERS' .or. satel(:3) == 'ENV') then
   irep=35000
else if (satel(:3) == 'CRY') then
	irep=369000
else
   stop 'Unknown satellite'
endif
call setearth(6378136.3d0,298.257d0)

! Open ODR file

open (20,file=odrname,status='replace',form='unformatted',recl=16,access='direct')
i=index(odrname,'ODR.')
read (odrname(i+4:i+6),*) iarc

! Process POE data records

600 format (6d22.16)
do
	read (*,600,iostat=ios) time,sec,grhran,xpole,ypole
	if (ios /= 0) exit
	read (*,600) xyz
	read (*,600) ecf
	read (*,*)

! Store date and start-time

	call rotate(1,-ypole*masrad,ecf,ecf)
	call rotate(2,-xpole*masrad,ecf,ecf)
	call xyzgeo(ecf,r,dlat,dlon,dhgt)
	dlat=dlat/rad
	dlon=dlon/rad
	if (dlon > 180d0) dlon=dlon-360
	krec=krec+1

	o%itime=nint(sec85(4,time*1d2+sec))
	o%lat=nint(dlat*1d7)
	o%lon=nint(dlon*1d7)
	o%hgt=nint(dhgt*1d3)
	write (20,rec=krec+2) o
enddo

! Reread the first and the last record.

read (20,rec=3) itime0
read (20,rec=krec+2) itime1
interval = nint((itime1-itime0)/dble(krec-1))

! Check if there was a leapsecond in the data period
! If so, determine when the leapsecond occurred, and correct all
! timetags after the event.
! This is because of the block structure of the G2T file.

if (mod(itime0,interval) /= mod(itime1,interval)) then
	mjd=itime1/86400+mjd85
	mjd=mdate(2,mdate(1,mjd)/100*100+01)
	leapsec=(mjd-mjd85)*86400
	do i=3,krec+2
		read (20,rec=i) o
		if (o%itime > leapsec .and. mod(o%itime,interval) == mod(itime0,interval)) then
			o%itime=o%itime-1
			write (20,rec=i) o
		endif
	enddo
endif

! Search for begin of precise arc.
! If skip = 0: start immediately
! If skip < 0: start immediately only when arc does not start at noon or midnight
! Else: start around abs(skip) days from beginning or arc

if (begin == 0d0) begin = itime0
if (skip == 0d0) then
	skip = begin
else if (modulo(itime0,43200) /= 0 .and. skip < 0d0) then
	skip = begin
else
	skip = begin + abs(skip)*86400
	call search(dble(itime0),interval,skip,krec)
endif
iskip = nint(skip)

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
real(eightbytereal) time0,skip,dlat,x
integer(fourbyteint) :: step,krec,irec,itime0,itime1,itime2,lat0,lat1,lat2

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

end program poe2odr
