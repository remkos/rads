program spx2odr

! Convert SP1, SP3 or SP3c file to Orbital Data Records (ODR).
!
! syntax: spx2odr odrfile < spxfile
!-
! $Log: spx2odr.f90,v $
! Revision 1.17  2020/11/10 09:52:46  rads
! - Add Sentinel-6
!
! Revision 1.16  2018/02/28 13:28:14  rads
! - Add Sentinel-3
!
! Revision 1.15  2017/11/22 16:10:21  rads
! - Allow for timing of orbit files different from integer seconds
!
! Revision 1.14  2017/03/21 12:56:40  rads
! - Update to support SP3d format
!
! Revision 1.13  2013/07/02 20:31:47  rads
! - Added support for SARAL
!
! Revision 1.12  2012/11/01 05:23:17  rads
! - Allow either ODR.xxx or arc_xxxx.odr
!
! Revision 1.11  2012/10/29 19:21:43  rads
! - Need to clip files at leap second when input jumps by one second
!
! Revision 1.10  2012/06/11 19:15:10  rads
! - Esthetics only
!
! Revision 1.9  2011/09/09 16:11:08  rads
! - Improve working with concatenated SP3 files
!
! Revision 1.8  2011/01/24 13:24:29  rads
! - Added CRYOSAT2 support
! - Allow TAI time system
! - Allow overruling of time system
!
! Revision 1.5  2010/12/13 21:05:50  rads
! - Added begin= option
!
! Revision 1.4  2010/11/12 01:16:49  rads
! - Replace ODR when it exists
!
! Revision 1.3  2010/06/14 19:16:45  rads
! - Application of leap second could be off by one record
!
! Revision 1.2  2010/05/29 00:37:01  rads
! - Handle orbits that do not run at :00 or :30 seconds
! - Change the way skip is dealt with
!
! Revision 1.1  2010/05/28 16:49:46  rads
! - Based on sp12odr.f by Eelco Doornbos and Remko Scharroo
! - Converted from Fortran 77 to Fortran 99
! - Added support for SP3 and SP3c
! - Convert all orbits to TOPEX ellipsoid, also from ERS and Envisat
!
! Created by Remko Scharroo - Altimetrics LLC
!***********************************************************************
use typesizes
character(80) :: odrname=' ',arg,altim
character(4) :: spec='xODR'
character(8) :: satel='ERS-1'
integer(fourbyteint) :: lstr,i,krec=0,itime0,itime1,iargc,irem,irep,iarc,iskip,mjd,ios,offset
real(eightbytereal) :: ecf(3),ecf0(3)=0d0,utc,t,skip=-1d0,dlon,dlat,dhgt,r,begin=0d0,sec85
integer(fourbyteint), parameter :: mjd85=46066
real(eightbytereal), parameter :: pi=4d0*atan(1d0), rad=pi/180d0
logical :: ltlend

! Variables for ODR format

type :: odr
	integer(fourbyteint) :: itime,lat,lon,hgt
endtype
type(odr) :: o

! Variables for SPx format

character(80) :: line
integer(fourbyteint) :: yyyy, month, day, hour, minute, mjdstart, num, numsat
real(eightbytereal) :: interval, fract, seconds
character(1) :: sp_version
character(3) :: orbtype
character(4) :: agency
character(5) :: dataused,coordsys
character(3) :: timesys = ''
logical :: skip_next = .false.

! Variables for TAI

integer(fourbyteint), parameter :: nutc=100
real(eightbytereal) :: tai_utc(nutc)
integer(fourbyteint) :: mjdutc(nutc)
real(eightbytereal), parameter :: tai_gpst=19d0
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
	else if (arg(:6) == 'begin=') then
		read (arg(7:),*) begin
		begin = sec85(0,begin)
	else if (arg(:4) == 'sys=') then
		timesys=arg(5:)
	else
		odrname=arg
	endif
enddo

! Write usage info on faulty input

if (odrname == ' ') then
	write (*,1301)
1301 format ('spx2odr: Convert SP1 or SP3 orbit file to ODR file'// &
'syntax: spx2odr [options] < spxfile odrfile'// &
'where [options] are:'/ &
'  sat=SATEL   : specify 8-character satellite name (def: GFO-1)'/ &
'  rem=REMID   : specify orbit ID'/ &
'  skip=DAY    : specify days to skip at start of orbit'/ &
'  begin=TIME  : begin arc at given time')
	stop
endif

! Check satellite

if (satel(:3) == 'GFO') then
	irep=17050
else if (satel(:3) == 'TOP' .or. satel(:3) == 'JAS' .or. satel(:7) == 'SNTNL-6') then
	irep=9920
else if (satel(:3) == 'ERS' .or. satel(:3) == 'ENV' .or. satel(:3) == 'SAR') then
	irep=35000
else if (satel(:3) == 'CRY') then
	irep=369000
else if (satel(:7) == 'SNTNL-3') then
	irep=27000
else if (satel(:3) == 'CHA') then
	irep=0
else
	stop 'Unknown satellite'
endif
call setearth(6378136.3d0,298.257d0)

! Open ODR file

open (20,file=odrname,status='replace',form='unformatted',recl=16,access='direct')
i=index(odrname,'ODR.')
if (i > 0) read (odrname(i+4:i+6),*) iarc
i=index(odrname,'.odr')
if (i > 0) read (odrname(i-4:i-1),*) iarc

! Read first SPx header record

read (*,'(a80)') line
if (line(1:3) == ' # ' .or. line(1:3) == '_#_') then	! SP1
	read (line,'(4x,i4,4i3,f11.7,f15.7,i6,f16.13,i7,a1,a4)') &
	yyyy,month,day,hour,minute,seconds,interval,mjdstart,fract,num,orbtype,agency
	sp_version='1'
	read (*,'(3x,i2)') numsat
	if (numsat > 1) stop 'Multiple satellites in SP1'
else if (line(1:1) == '#' .and. line(2:2) >= 'a' .and. line(2:2) <= 'd') then	! SP3
	read (line,'(3x,i4,4i3,f12.8,i8,1x,a5,1x,a4,1x,a3,1x,a4)') &
	yyyy,month,day,hour,minute,seconds,num,dataused,coordsys,orbtype,agency
	read (*,'(23x,f15.8,i6,f16.13)') interval,mjdstart,fract
	read (*,'(4x,i2)') numsat
	if (numsat > 1) stop 'Multiple satellites in SP3'
else
	stop 'Invalid SP1/SP3 header line 1'
endif

! Check if we have an offset from integer seconds

offset = nint((seconds - nint(seconds)) * 1d6)
if (offset /= 0) then
	irem = offset - 1000000
	write (*,"('Time offset detected: ',f9.6)") offset * 1d-6
endif

! Read SPx epoch and position/velocity records

do
	read (*,'(a80)',iostat=ios) line
	if (ios /= 0) exit

	! First look for epoch record

	if (line(1:3) == 'EOF') then
		cycle
	else if (line(1:2) == '%c' .and. line(10:12) /= 'ccc') then
		! Set time system if not undefined
		timesys = line(10:12)
		cycle
	else if (line(1:3) == ' # ' .or. line(1:3) == '_#_') then
		cycle
	else if (line(1:1) == '+' .or. line(1:1) == '/' .or. line(1:1) == '%' .or. line(1:1) == '#') then	! SP3 headers
		! If headers are interspersed, skip the next state vector
		if (krec > 0) skip_next = .true.
		cycle
	else if (line(1:3) == ' + ' .or. line(1:3) == '_+_') then
		! If headers are interspersed, skip the next state vector
		if (krec > 0) skip_next = .true.
		cycle
	else if (sp_version == '1' .and. (line(1:3) == ' * ' .or. line(1:3) == '_*_')) then
		read (line,'(4x,i4,4i3,f11.4)') yyyy,month,day,hour,minute,seconds
	else if (line(1:2) == '* ') then
		read (line,'(3x,i4,4i3,f11.8)') yyyy,month,day,hour,minute,seconds
	else
		write (*,*) line
		stop 'Epoch record expected'
	endif

	! Then read state vector

	read (*,'(a80)',iostat=ios) line
	if (ios /= 0) exit
	if (line(1:1) == 'P') then
		read (line,'(4x,3f14.6)') ecf
		read (*,*)  ! Skip velocity
	else if (line(1:4) == ' SV ') then
		read (line,'(6x,3f13.6)') ecf
	else if (line(1:3) == ' SV') then
		read (line,'(5x,3f13.6)') ecf
	else
		stop 'State vector record expected'
	endif
	ecf = ecf*1d3   ! km to m

	! If state vector is to be skipped (after any but first header) ...
	if (.not.skip_next) then
	else if (dot_product(ecf-ecf0,ecf-ecf0) > 1d6) then
		! Skipped state vectors should be the same as the previous. If not, then there is a leapsecond and we should quit
		write (*,"('Leapsecond, file clipped at ',i4.4,'-',i2.2,'-',i2.2)") yyyy,month,day
		exit
	else
		skip_next = .false.
		cycle
	endif
	ecf0 = ecf

! Convert GPS time to UTC
! UTC = GPST + (TAI - GPST) - (TAI - UTC)
! or UTC = TAI - (TAI - UTC)

	call ymd2mjd (yyyy,month,day,mjd)
	t = (mjd-mjd85)*86400d0+hour*3600d0+minute*60d0+seconds
	if (timesys == 'GPS') then
		utc = t + tai_gpst - tai_utc(iutc)
		do while (utc/86400d0 >= mjdutc(iutc+1) .and. modulo(utc,86400d0) > 0.5d0 .and. iutc < nutc)
			iutc = iutc + 1
		enddo
		utc = t + tai_gpst - tai_utc(iutc)
	else if (timesys == 'TAI') then
		utc = t - tai_utc(iutc)
		do while (utc/86400d0 >= mjdutc(iutc+1) .and. modulo(utc,86400d0) > 0.5d0 .and. iutc < nutc)
			iutc = iutc + 1
		enddo
		utc = t - tai_utc(iutc)
	else ! UTC
		utc = t
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

end program spx2odr
