      program asc2gbx
*
* Converts ASCII records to GEODYN I Binary Extended for altimetry
*
* ascii records should contain:
* mjd    integer: modified julian day
* sec    real:    seconds of the day
* dlat   real:    latitude (deg)
* dlon   real:    longitude (deg)
* h      real:    altimeter range measurement (m)
* sigmah real:    altimeter range error (m)
*
* Remko Scharroo - DUT/SSR&T - 12 Jan 1994.
*
      implicit none
      real*8    dayfrac,rad
      real*8    sec,dlat,dlon,h,dhellips,sig/0d0/,dlon0
      integer*4 gbx4(27)/27*0/,b(0:31),mjd,itide,isst,i,n,ikey,isat,
     |		iargc,test/0/
      integer*2 gbx2(54)
      real*4    sigmah
      character*80 filenm,arg
      equivalence (gbx2,gbx4),(gbx4(5),mjd),(gbx4(12),sigmah)
      parameter (itide=3,isst=1)
      real*8	ae_grs,f_grs,ae_new,f_new,amp/1d-1/
      parameter (ae_grs=6378137.0d0,f_grs=1/298.2570d0)
      parameter (ae_new=6378136.3d0,f_new=1/298.2564d0)
      logical ellips/.true./

      filenm=' '
      sig=0d0
      do i=1,iargc()
	 call getarg(i,arg)
	 if (arg(:4).eq.'sig=') then
	    read (arg(5:),*) sig
	 else if (arg(:5).eq.'test=') then
	    read (arg(6:),*,iostat=n) test,amp
	 else if (arg(:2).eq.'-e') then
	    ellips=.true.
	 else
	    filenm=arg
	 endif
      enddo
      if (filenm.eq.' ') goto 1300

* Initialize

      gbx2(3)=49	! measurement type (49=xover, 43=alt)
      gbx2(4)=23	! time system
      gbx4(3)=0		! station number
      gbx4(13)=0	! net instrument correction
      gbx4(15)=0	! meteorological word
      gbx4(16)=0	! geoid height
      gbx4(17)=0	! ocean dynamic correction
      gbx4(18)=0	! surface elevation
      gbx4(19)=0	! spacecraft revolution number
      gbx4(20)=0	! mean sea surface elevation
      gbx2(43)=0	! H1/3
      gbx2(44)=0	! AGC
      gbx2(45)=0	! wind speed
      gbx2(47)=0	! dry tropo
      gbx2(48)=0	! wet tropo
      gbx2(50)=0	! iono
      gbx2(51)=0	! barotropic dynamic sea surefcae correction
      gbx2(52)=0	! solid earth tides
      gbx2(53)=0	! Schwiderski ocean tides
      gbx2(54)=0	! other ocean tide

      rad=atan(1d0)/45

      open (20,file=filenm,status='new',form='unformatted')
      n=0
      do i=0,31
         b(i)=2**i
      enddo
      gbx4(4)=3*b(25)+b(24)+itide*b(22)+b(21)+b(20)+b(19)
      gbx2(46)=b(14)+b(13)+isst*b(12)+b(11)+b(10)+b(8)+
     |         b(7)+b(6)+b(5)+b(3)+b(2)+b(1)+b(0)

* Read input records

  100 read (5,*,iostat=i,end=900) mjd,sec,dlat,dlon,h,ikey,sigmah,isat

      if (ikey.eq.0) gbx2(3)=43	! If key=0 then direct altimeter meas.

      if (isat.ne.0) then
	 gbx4(1)=isat		! copy satellite ID from record
      else if (h.lt.1d6) then
         gbx4(1)=9105001	! ERS-1
      else
         gbx4(1)=9205201	! TOPEX/Poseidon
      endif
      dayfrac=sec/86400		! day fraction
      gbx4(10)=nint(dlat*1d6)	! latitude
      gbx4(11)=nint(dlon*1d6)	! longitude
      gbx4(21)=ikey		! xover key
      sigmah=sqrt(sigmah**2+sig**2)
      if (sig.lt.0 .or. sigmah.lt.0) sigmah=-sigmah

* Convert height from Geosat (6378137.0,298.257) to
* TOPEX (6378136.3,298.2564) ellipse

      if (ellips) h=h+dhellips(3,dlat)

* For test purpuses: add one of the error models

      if (test.eq.1) then
         h=h+amp*sin(dlat*rad)
      else if (test.eq.2) then
         h=h+amp*cos(dlat*rad)*cos(dlon*rad)
      else if (test.eq.3) then
         h=h+amp*cos(dlat*rad)*sin(dlon*rad)
      else if (test.eq.4) then
	 dlon0=(sec-81d3)/240d0
         h=h+amp*cos(dlat*rad)*cos((dlon+dlon0)*rad)
      endif

      n=n+1
*     write (20) gbx4
      write (20) (gbx4(i),i=1,5),dayfrac,h,(gbx4(i),i=10,27)
*     write (20) gbx4(1),(gbx2(i),i=3,4),
*    |           (gbx4(i),i=3,5),dayfrac,h,
*    |		 (gbx4(i),i=10,20),(gbx2(i),i=41,54)
      goto 100

 1300 write (0,1301)
 1301 format ('usage: asc2gbx [options] gbxfile'//
     |'Options are:'/
     |'sig=sigma  : add sigma (in RSS sense) to value on records'//
     |'Reads ascii records (t,lat,lon,h,key,sigma,satid) from stdin'/
     |'and writes GBX file')
      goto 999

  900 write (0,910) n
  910 format (i6,' records converted asc -> gbx')
  999 end
