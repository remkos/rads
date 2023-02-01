**STATINFO -- Return information for given laser station
*+
      FUNCTION STATINFO (TIME, STATION, SITENAME, OCCODE, CODE,
     |			 NOISE, WAVELENGTH, PLATE, ECCENTRICITY)
      REAL*8	TIME, ECCENTRICITY(3)
      CHARACTER SITENAME*32, OCCODE*6, CODE*4
      INTEGER*4 STATION, NOISE, WAVELENGTH, PLATE, STATINFO

* For given STATION and TIME, this routine returns the correct WAVELENGTH,
* SITENAME, occupation CODE, and system NOISE. When TIME is >= 1d20 then
* TIME is not used as a selection criterion. The most recent information for
* STATION is returned.
*
* The routine initialises upon first the first call. It requires the access
* to file ${ALTIM}/data/tables/system.data.
*
* The routine uses a fast query mechanism to get the station information
* with only a little expense of memory.
*
* Arguments:
*  TIME        (input) : Time in MJD
*  STATION     (input) : Station ID
*  SITENAME   (output) : Name of the laser station site
*  OCCODE     (output) : 6-letter occupation code
*  CODE       (output) : 4-letter IGS code
*  NOISE      (output) : Laser ranging noise level in centimetres
*  WAVELENGTH (output) : Wavelength of laser station in nanometres
*  PLATE      (output) : Tectonic plate ID
*  ECCENTRICITY  (out) : North, East, and Up component of eccentricity
*                        vector in metres.
*
* Exit code:
*  STATINFO   (output) : 0=no error, 1=no entry found
*-
*  5-Aug-1998 - Remko Scharroo
*  9-Aug-1998 - Increase number of returned variables. Accellerated searching.
* 13-Apr-2001 - Bug fix: include SAVE plat,ecce
*  1-May-2001 - Added IGS code as output
* 28-Feb-2004 - Multiple warnings suppressed. nstat saved (finally)
*-----------------------------------------------------------------------
      integer*4 mstat,nstat,freeunit,mdate,ioerr
      parameter (mstat=1500)
      integer*4 pnt1(10000),pnt2(10000),warn(10000)
      integer*4 stat(mstat),t0(mstat),t1(mstat),wave(mstat),nois(mstat),
     |		plat(mstat)
      real*4	ecce(3,mstat)
      character name(mstat)*32,cod1(mstat)*6,cod2(mstat)*4
      integer*4 i,j,ios,unit,itime
      character filename*80
      save nstat,stat,t0,t1,wave,nois,plat,ecce,name,cod1,cod2,
     |		pnt1,pnt2,warn
      data nstat/0/

* Upon first call: read the system.data file and keep all info regarding
* site occupation, system noise, eccentricity vector, etc.
* Also initialise an array with the location of the first entry for any
* given station ID (PNT1)

      if (nstat.eq.0) then
         do i=1,10000
	    pnt1(i)=0
	    warn(i)=0
	 enddo
         unit=freeunit()
         filename = ' '
         call checkenv('SYSTEMDATA',filename,i)
         if(filename(1:1).eq.' ') then
           filename='/user/altim'
           call checkenv('ALTIM',filename,i)
           filename(i+1:)='/data/tables/system.data'
         end if
         open (unit,file=filename,status='old',iostat=ioerr)
         if (ioerr .ne. 0) THEN
           write(0,*) "Error in statinfo: Can not open ", filename
           stop
         end if
         rewind (unit)
10       nstat=nstat+1
11	 if (nstat.gt.mstat)
     |		call fin('statinfo: too many entries in system.data')
         read (unit,20,iostat=ios,end=100) stat(nstat),cod1(nstat),
     |		(ecce(i,nstat),i=1,3),t0(nstat),t1(nstat),wave(nstat),
     |		plat(nstat),name(nstat),nois(nstat),cod2(nstat)
20       format (i4,6x,a6,3f8.3,2i8,i5,6x,i4,7x,a32,i2,t131,a4)
	 if (ios.ne.0) goto 11
	 t0(nstat)=mdate(2,t0(nstat))
	 t1(nstat)=mdate(2,t1(nstat))
	 i=stat(nstat)
	 if (pnt1(i).eq.0) pnt1(i)=nstat
	 pnt1(i+1)=nstat+1
         goto 10

100      close (unit)
         nstat=nstat-1
	 do i=1,10000
	    pnt2(i)=pnt1(i)
	 enddo
      endif

* Once initialised, we try to jump directly to the correct line.
* PNT2 is the location of the last visited line for a given station ID
* PNT1(i) points to the first line for a given station ID (i) and
* PNT1(i+1)-1 is the last line.
* When time >= 1d20, then give the last line.

      if (station.lt.0 .or. station.gt.9999) then
	 statinfo=1
	 return
      endif
 
      if (time.ge.1d20) then
	 i=pnt1(station+1)-1
	 if (stat(i).eq.station) goto 900
      endif

      itime=int(time)
      do i=pnt2(station),pnt1(station+1)-1
         if (stat(i).eq.station .and.
     |     itime.ge.t0(i) .and. itime.le.t1(i)) goto 900
      enddo

      do i=pnt1(station),pnt2(station)-1
         if (stat(i).eq.station .and.
     |     itime.ge.t0(i) .and. itime.le.t1(i)) goto 900
      enddo

* If we get here, the correct station info was not found.

      wavelength=600
      noise=99
      sitename='(Unknown)'
      occode='******'
      code='****'
      plate=0
      eccentricity(1)=0
      eccentricity(2)=0
      eccentricity(3)=0
      statinfo=1

      if (warn(station).eq.0) then
	 write (0,600) station,time
      else if (mod(int(warn(station)),10).eq.0) then
	 write (0,601) station,warn(station)
      endif
      warn(station)=warn(station)+1
600   format ('STATINFO: WARNING:',
     |   ' no entry could be found for station',i5,
     |   ' for epoch',f10.3,
     |   '. Default values are assumed.')
601   format ('STATINFO: WARNING for station',i5,
     |   ' repeated',i5,' times.')
      return

900   pnt2(station)=i
      wavelength=wave(i)
      noise=nois(i)
      sitename=name(i)
      occode=cod1(i)
      code=cod2(i)
      plate=plat(i)
      do j=1,3
         eccentricity(j)=ecce(j,i)
      enddo
      statinfo=0

      return
      end
