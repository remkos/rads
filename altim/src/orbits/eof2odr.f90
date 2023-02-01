program eof2odr

! Convert Earth Explorer Orbit File (EOF) to Orbital Data Records (ODR).
!
! syntax: eof2odr [ options ] odr_file < eof_file
!-
! $Log: eof2odr.f90,v $
! Revision 1.2  2020/11/27 16:30:20  rads
! - Added support for Sentinel-6
!
! Revision 1.1  2016/05/22 13:20:44  rads
! eof2odr first version, based on spx2odr
!
!
! Created by Remko Scharroo - EUMETSAT
!***********************************************************************
use typesizes
character(80) :: odrname=' ',arg
character(4) :: spec='xODR'
character(8) :: satel=''
integer(fourbyteint) :: i,n,krec=0,itime0,itime1,iargc, &
		irem=0,irep=0,iarc=0,iskip=0
real(eightbytereal) :: ecf(3),utc,utc0,utc1,skip=-1d0,dlon,dlat,dhgt,r
real(eightbytereal), parameter :: pi=4d0*atan(1d0), rad=pi/180d0, mjd85=46066d0
logical :: ltlend

! Variables for ODR format

type :: odr
	integer(fourbyteint) :: itime,lat,lon,hgt
endtype
type(odr) :: o

! Variables for EOF format

character(132) :: line
real(eightbytereal) :: interval = 60d0

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
1301 format ('eof2odr: Convert Earth Explorer Orbit File (EOF) to ODR file'// &
'syntax: eof2odr [options] < eof_file odr_file'// &
'where [options] are:'/ &
'  sat=SATEL   : specify satellite'/ &
'  rem=REMID   : specify orbit ID'/ &
'  skip=DAY    : specify days to skip at start of orbit')
	stop
endif

! Open ODR file

open (20,file=odrname,status='replace',form='unformatted',recl=16,access='direct')
i = index(odrname,'ODR.')
if (i > 0) read (odrname(i+4:i+6),*) iarc
i = index(odrname,'.odr')
if (i > 0) read (odrname(i-4:i-1),*) iarc
call setearth (6378136.3d0,298.257d0)

! Read EOF header records

call lees ('Mission', line)
if (index(line,'Sentinel-3A') > 0) then
	satel='SNTNL-3A'
	irep=27000
else if (index(line,'Sentinel-3B') > 0) then
	satel='SNTNL-3B'
	irep=27000
else if (index(line,'Sentinel-6A') > 0) then
	satel='SNTNL-6A'
	irep=9920
else if (index(line,'Sentinel-6B') > 0) then
	satel='SNTNL-6B'
	irep=9920
endif

call lees ('Validity_Start', line)
call time_convert (line, utc0)

call lees ('Validity_Stop', line)
call time_convert (line, utc1)

call lees ('List_of_OSVs', line)
read (line(8:len_trim(line)-1), *) n
interval = (utc1-utc0)/(n-1)

! Read OSV epoch and position/velocity records

do i = 1,n
	call lees ('UTC', line)
	call time_convert (line, utc)
	call lees ('X', line)
	read (line, *) ecf(1)
	call lees ('Y', line)
	read (line, *) ecf(2)
	call lees ('Z', line)
	read (line, *) ecf(3)

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

subroutine lees (tag, result)
character(len=*), intent(in) :: tag
character(len=*), intent(out) :: result
character(len=160) :: line
integer :: i, j

do
	read (*,'(a)') line
	i = index(line,'<'//tag)
	if (i > 0) exit
enddo
j = index(line,'</')
if (j == 0) then
	result = line(i+len_trim(tag)+2:len_trim(line)-1)
else
	i = index(line,'>')
	result = line(i+1:j-1)
endif
end subroutine lees

subroutine time_convert (line, utc)
character(len=*), intent(in) :: line
real(eightbytereal), intent(out) :: utc
integer(fourbyteint) :: yyyy,month,day,hour,minute,mjd
real(eightbytereal) :: seconds
read (line,'(4x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f9.6)') yyyy,month,day,hour,minute,seconds
call ymd2mjd (yyyy,month,day,mjd)
utc = (mjd-mjd85)*86400d0+hour*3600d0+minute*60d0+seconds
end subroutine time_convert

subroutine search(time0,step,skip,krec)
real(eightbytereal) :: time0,step,skip,dlat,x
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

end program eof2odr
