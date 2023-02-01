**GETORB -- Get orbital position of satellite from ODR files
*+
      FUNCTION GETORB (TIME, LAT, LON, ORBIT, PATH, VERBOSE)
      INTEGER*4	GETORB
      REAL*8	TIME, LAT, LON, ORBIT
      CHARACTER	PATH*(*)
      LOGICAL	VERBOSE
*
* This function reads the Orbital Data Records (ODR) in a directory or a
* file indicated by the character PATH and interpolates the position of
* the sub-satellite point (LAT and LON) and the orbital altitude (ORBIT)
* above the TOPEX reference ellipsoid ( ae = 6378137.3 meter; f = 1/298.257 )
* at the requested epoch TIME (in UTC seconds since 1985).
*
* If PATH starts with + the rest of the character string is assumed to be
* the pathname of a regular file. Otherwise, it is assumed to be a
* directory name.
*
* (1) PATH = directory
*
* If this function is called for the first time, it will scan the directory for
* existing ODR files (with the regular numbering ODR.iii), and store
* the begin and end epoch of the precise part of each file. This scanning
* is only repeated when a change of directory is made. The scanning is
* accelerated if a file 'arclist' is available in the directory with the
* same information.
*
* This routine will forthwith load the ODR file that contains the orbital
* positions around the requested epoch. A new file is only loaded when a
* next call of GETORB is made for an epoch TIME outside the limits of the
* loaded file, or when the directory path is changed.
* Therefore, GETORB works fastest when subsequent calls are made for epochs
* which are time sorted (ascending or descending), or at least sorted by
* arc number.
*
* This routine also supports orbital arcs including a prediction. A
* position will be interpolated up to the end of the last available ODR
* in the given directory.
*
* (2) PATH = +file
*
* A single ODR file is loaded at the first call with a new string PATH.
* The routine will not scan for alternative ODR files. Interpolation is
* restricted to this single ODR file indicated by the string following
* the plus sign.
*
* In both cases an error return code (GETORB > 0) is given it the requested
* epoch is outside the limits of the loaded file or the available files in
* the given directory. A warning return code (GETORB < 0) is given 
* when the file is interpolated outside the recommended precise part of
* the arc.
*
* Finally the position and altitude is interpolated. This interpolation is
* done in pseudo-Cartesian coordinates: latitude, longitude, and altitude
* are treated as polar coordinates and are changed to pseudo-Cartesian,
* and are converted back after interpolation of these coordinates.
*
* Arguments:
*  TIME    (input): Epoch at which the position of the satellite is requested.
*                   This epoch has to be given in UTC seconds since 
*                   1.0 January 1985.
*  LAT    (output): Geodetic latitude of the sub-satellite point at the
*                   requested epoch (in degrees).
*  LON    (output): Longitude of the sub-satellite point at the requested
*                   epoch (in degrees, interval -180 to 180 degrees).
*  ORBIT  (output): Orbital altitude above the TOPEX ellipsoid (in meters),
*                   interpolated at the requested epoch.
*  PATH    (input): (1) Name of the directory in which the ODRs are stored,
*                       e.g. '/home/ers1/ODR/DGM-E04'
*                   (2) Name of the ODR file, preceeded by a +
*  VERBOSE (input): Gives some processing information if .TRUE.
*  GETORB (output):     0 = no error/warning
*                   Fatal error codes:
*                       1 = error opening ODR file,
*                       2 = no ODR for requested epoch, 3 = too many records,
*                   Warning codes (non fatal):
*                      -1 = time within ODR, but outside precise part.
*-
* - Where to get orbit files ?
* ODR orbit files and arclists for ERS-1 and ERS-2 are to be obtained from
* the anonymous-ftp server of NOAA Geophysics Lab: harpo.grdl.noaa.gov,
* directories pub/delft/ODR.ERS-?/<model> where <model> relates to the
* gravity field model and/or additional tracking data used in the orbit
* computation.
*
* - Supported platforms
* This subroutine automatically recognises Little and Big Endian representation
* of integers, and hence is able to convert the input files to the proper
* representation.
*
* For more information, e-mail to Remko Scharroo, DUT/DEOS.
* remko.scharroo@lr.tudelft.nl
*
* This routine requires: INTER8, INTAB8, FREEUNIT, I4SWAP,
* POLCAR, CARPOL.
*-
* $Log: getorb.f,v $
* Revision 1.10  2017/11/22 16:10:21  rads
* - Allow for timing of orbit files different from integer seconds
*
* Revision 1.9  2015/07/10 14:29:59  rads
* Refer to TOPEX, not GRS80 ellipsoid.
*
* Revision 1.8  2014/01/07 16:30:00  rads
* - Deal better with gaps between orbit files.
* - Properly indicate that orbit is not valid during gap.
*
* Revision 1.7  2010/11/15 18:10:05  rads
* - Expanded to 4-digit ODR file numbering (arc_????.odr)
*
* Revision 1.6  2010/06/14 21:08:31  rads
* - Send messages to STDOUT instead of STDERR
*
* Revision 1.5  2007/03/26 20:25:54  rads
* - Check outer time limits before inner time limits
*
* Revision 1.4  2006/07/28 22:03:53  rads
* - Removed obsolete Fortran constructions
*
* Revision 1.3  2006/01/10 14:19:15  rads
* - Removed FASTIO; reverted to FORTRAN I/O
*
* Revision 1.2  2004/11/22 14:30:26  remko
* - Allow arc numbers higher than 999
*
*  9-Aug-1999 - Removed closef when odrinfo fails
* 12-Nov-1996 - New ODR format introduced. Less code.
* 22-Feb-1996 - Handle leap seconds properly.
* 18-May-1995 - Multi-platform version.
* 11-Jan-1994 - Scan arclist to speed up arc selection
* 18-Jun-1993 - GETORBFILE included
* 18-Nov-1992 - Revised version.
*------------------------------------------------------------------------------
      real*8   	rev
      integer*4	maxnrarcs
      parameter	(maxnrarcs=9999)
      integer*4	istart,iend,iskip(maxnrarcs,2),iarc,iskip0,iskip1,
     |	itime0,itime1,itstep,getorbfile,odrinfo,mdate,lnblnk,
     |	i,j,l,jarc,irep,nrec,remid,maxarc,minarc,nrev,unit,freeunit,
     |  yymmdd0,hh0,mm0,yymmdd1,hh1,mm1,yymmdd2,hh2,mm2,ss2,mode/0/
      logical*4	exist(0:maxnrarcs)
      character	satel*8,odrnm*80,old_path*80,line*80,arc*4

      include	'math.inc'

      common	/cgetorb/ rev,nrev,satel
      save
      data	old_path/' '/

* Directory or file

      if (path(1:1).eq.'+') then
	 getorb=getorbfile(time,lat,lon,orbit,path(2:),verbose)
	 return
      endif

* What to do if new path is selected

      if (path.ne.old_path) then

	 old_path=path
	 l=lnblnk(path)
	 maxarc=0
	 minarc=maxnrarcs

* Check if arclist exists.

	 odrnm=path
	 odrnm(l+1:)='/arclist'
	 inquire (file=odrnm,exist=exist(0))
	 if (exist(0)) then

* Arclist exists. Scan arclist.

	    if (verbose) write (*,608) odrnm(:l+8)
	    unit=freeunit()
	    open (unit,file=odrnm,status='old')
	    rewind(unit)
9	    read (unit,550,end=11) line
	    if (line(1:4).ne.'Arc#') goto 9
10	    read (unit,620,end=11) arc,yymmdd0,hh0,mm0,
     |		yymmdd1,hh1,mm1,yymmdd2,hh2,mm2,ss2
	    if (arc(4:4) .eq. ' ') then
	    	read (arc,'(z1,i2)') i,j
		iarc = i*100 + j
		mode = 0
	    else
	    	read (arc,'(i4)') iarc
		mode = 1
	    endif
	    exist(iarc)=.true.
	    itime0=(mdate(2,yymmdd0)-46066)*86400+hh0*3600+mm0*60
	    itime1=(mdate(2,yymmdd1)-46066)*86400+hh1*3600+mm1*60
	    istart=(mdate(2,yymmdd2)-46066)*86400+hh2*3600+mm2*60+ss2
	    iskip(iarc,1)=istart
	    iskip(iarc,2)=itime1
	    if(iarc.gt.1 .and. istart.lt.iskip(iarc-1,2))
     |	      iskip(iarc-1,2)=istart
	    minarc=min(minarc,iarc)
	    maxarc=max(maxarc,iarc)

* The first skip may be set to the beginning of the first arc

	    if (iarc.eq.minarc) iskip(iarc,1)=itime0
	    goto 10
11	    close (unit)
	    if (verbose) write (*,611)
	 else

* No arclist. Scan directory.

	    if (verbose) write (*,609) odrnm(:l)
	    do iarc=1,maxnrarcs-1
	       write (odrnm,1110) path(:l),iarc/100,mod(iarc,100)

* Check for existance of ODR. If headers are garbled, ignore file.

	       inquire (file=odrnm,exist=exist(iarc))
	       if (.not.exist(iarc)) goto 5
	       i=odrinfo(unit,odrnm,satel,irep,jarc,remid,
     |	   nrec,itime0,itime1,itstep,istart,iend,rev)
	       if (i.lt.0 .or. nrec.eq.0) then
		  exist(iarc)=.false.
	       else
		  close(unit)
		  iskip(iarc,1)=istart
		  iskip(iarc,2)=itime1
	          if(iarc.gt.1.and.istart.lt.iskip(iarc-1,2))
     |		    iskip(iarc-1,2)=istart
		  minarc=min(minarc,iarc)
		  maxarc=max(maxarc,iarc)

* The first skip may be set to the beginning of the first arc

	         if (iarc.eq.minarc) iskip(iarc,1)=itime0
	       endif
5	       continue
	    enddo
	    if (verbose) write (*,611)
	 endif

* iskip0 and iskip1 are current begin and end of precise arc.

	 iskip0=0
	 iskip1=0
      endif

* If epoch not between 'skips':
* - Check next skip
* - Load new orbit

      if (time.lt.iskip0 .or. time.gt.iskip1) then
	 do iarc=minarc,maxarc
	    if (time.ge.iskip(iarc,1) .and. time.le.iskip(iarc,2) .and.
     |	  exist(iarc)) goto 15
	 enddo
	 getorb=2
	 return

* Proper arc found. Continue with this one.

15	 continue
	 iskip0=iskip(iarc,1)
	 iskip1=iskip(iarc,2)
	 if (mode .eq. 0) then
	 	write (odrnm,1110) path(:l),iarc/100,mod(iarc,100)
      	 else
	 	write (odrnm,1111) path(:l),iarc
	 endif
      endif

* Call GETORBFILE

      getorb=getorbfile(time,lat,lon,orbit,odrnm,verbose)
      return

  550 format (a)
  608 format ('(getorb: reading arclist ',a,$)
  609 format ('(getorb: scanning directory ',a,$)
  611 format (')')
  620 format (a4,1x,i6,1x,i2,1x,i2,3x,i6,1x,i2,1x,i2,
     |          31x,i6,1x,i2,1x,i2,1x,i2)
 1110 format (a,'/ODR.',z1.1,i2.2)
 1111 format (a,'/arc_',i4.4,'.odr')

      end

      function getorbfile(time,lat,lon,orbit,file,verbose)

      character	old_file*80,file*(*),satel*8
      integer*4 maxnrrecs,mjd85
      parameter (maxnrrecs=60000,mjd85=46066)
      real*8	ovec(3,maxnrrecs),vec(3),time0,time1,tleap,
     |	rev,time,lat,lon,orbit,murad
      integer*4 getorbfile,itime,itime0,itime1,leap,type,
     |	ilat,ilon,iorbit,odr(4),l,irec,odrinfo,lnblnk,
     |	irep,iarc,remid,nrec,itstep,istart,iend,unit,mdate
      logical*4 verbose,swap,ltlend
      equivalence (odr(1),itime),(odr(2),ilat),
     |		(odr(3),ilon),(odr(4),iorbit)

      include 'math.inc'
      save
      data	old_file/' '/

      if (file.ne.old_file) then

* Open next orbit data file

	 old_file=file
	 l=lnblnk(file)
	 if (verbose) write (*,610) file(1:l)

* Open ODR file and get information from header.
* When old format, LAT and LON are in 10^-6 degrees. In new format
* they are in 10^-7 degrees.

	 type=odrinfo(unit,file,satel,irep,iarc,remid,
     |	nrec,itime0,itime1,itstep,istart,iend,rev)
	 if (type.lt.0) then
	    getorbfile=1
	    return
	 else if (nrec.gt.maxnrrecs) then
	    close(unit)
	    write (*,550) 'getorb: ODR file too big'
	    getorbfile=1
	    return
	 endif
         if (type.eq.0) then
            murad=rad/1d6
	 else
	    murad=rad/1d7
	 endif
	 swap=ltlend()

* Read all records and convert LAT,LON,ORBIT to pseudo-X,Y,Z

	 do irec=1,nrec
	    read (unit,rec=irec+2) odr
	    if (swap) call i4swap(4,odr)
	    call polcar(ilat*murad,ilon*murad,
     |		iorbit/1d3,ovec(1,irec))
	 enddo

* If leapsecond is in the arc, determine when it occurred and fiddle
* around a bit with the time interval.

	 if (mod(itime1-itime0,itstep).ne.0) then
	    leap=mdate(1,itime1/86400+mjd85)/100*100+01
	    if (verbose) write (*,612) leap
	    tleap=(mdate(2,leap)-mjd85)*86400
	    itime1=itime1+1
	 else
	    tleap=1d30
	 endif
	 time0=itime0
	 time1=itime1

* In case of offset in time (marked by remid < 0) then add the offset to time0 and time1

         if (remid.lt.0) then
	     time0 = time0 + (remid * 1d-6 + 1d0)
	     time1 = time1 + (remid * 1d-6 + 1d0)
	 endif

* Close the ODR file

	 if (verbose) write (*,611)
	 close(unit)
      endif

* Check time limits

      if (time.lt.itime0 .or. time.gt.itime1) then
	 getorbfile=2
	 return
      else if (time.ge.istart .and. time.le.iend) then
         getorbfile=0
      else
	 getorbfile=-1
      endif

* Interpolate to requested date, add one second if after leapsecond

      if (time.lt.tleap) then
         call intab8(3,nrec,ovec,time0,time1,time,vec)
      else
         call intab8(3,nrec,ovec,time0,time1,time+1d0,vec)
      endif

* Convert back to polar

      call carpol(vec,lat,lon,orbit)
      lat=lat/rad
      lon=lon/rad

      return

  550 format (a)
  610 format ('(getorb: orbit file ',a,$)
  611 format (')')
  612 format (' ... leapsecond: ',i6.6,$)

      end
