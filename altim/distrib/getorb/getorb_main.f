**GETORB -- An interface to the GETORB subroutine
*+
      program getorb_main

* This program functions as an interface to the GETORB subroutine.
* For a short syntax description run getorb without arguments.
*-
*    Dec-1996 - Created by Remko Scharroo - DEOS
* 25-Feb-1997 - XYZ added in output
*  9-Nov-1999 - Y2K compliancy
* 26-Nov-1999 - New manual
*  7-Dec-2000 - sysdep.inc removed
*-----------------------------------------------------------------------
      implicit  none
      real*8    time1,time2,step,time,lat,lon,orbit,sec85,xyz(3),
     |		r,rad,lon0/-180d0/,lon1/180d0/,lat0/-90d0/,lat1/90d0/
      integer*4 i,iarg,iargc,getorb
      character*80 path/' '/,arg
      logical*4 datearg,skip/.false./
      external  getorb

1310  format ('getorb -- a program to interpolate the ODR files'/,/
     |'usage: getorb [t=|ymd=|mjd=|doy=|sec=]t0[,t1,dt]',
     |' [lon=lon0,lon1] [lat=lat0,lat1] [+]path'/,/
     |'where:'/,/
     |'  t0      Epoch for the orbit interpolation in various formats:'/
     |'          ymd=t0: [YY]YYMMDD.DDD or [YY]YYMMDDHHMMSS.SSS'/
     |'          mjd=t0: Modified Julian Dates (MJD.DDD)'/
     |'          doy=t0: Year and day-of-year ([YY]YYDDD.DDD)'/
     |'          sec=t0: UTC seconds since 1 Jan 1985'/
     |'            t=t0: best guess between ymd= and mjd='/
     |'              t0: best guess between ymd= and mjd='/,/
     |'t0,t1,dt  Interpolation for a time range (t0 to t1) requested'/
     |'          (formats as above) with step size dt (in seconds)'/,/
     |'lon0,lon1 longitude range (deg)'/
     |'lat0,lat1 latitude range (deg)'/,/
     |'  path    path name of directory in which ODR files are stored'/
     |' +path    path name of a single ODR file to be used'/,/
     |'Examples:'/
     |'   getorb mjd=50001.500 /home/ers2/ODR'/
     |'   getorb ymd=991231235959,000101000159,60 .')

* Initialize

      rad=atan(1d0)/45
      time1=1d40
      time2=1d40
      step=1d40

* Read the time argument

      do iarg=1,iargc()
         call getarg(iarg,arg)
         if (datearg(arg,time1,time2,step)) then
         else if (arg(:4).eq.'lat=') then
            read (arg(5:),*,iostat=i) lat0,lat1
         else if (arg(:4).eq.'lon=') then
            read (arg(5:),*,iostat=i) lon0,lon1
         else if (time1.gt.1d30) then
            read (arg,*,iostat=i) time1,time2,step
            time1=sec85(0,time1)
         time2=sec85(0,time2)
         else
      	    path=arg
         endif
      enddo

      if (time1.gt.1d40 .or. path.eq.' ') then
c         write (*,*) path,time1
         write (*,1310)
	 goto 9999
      endif

      if (time2.gt.1d30) time2=time1
      if (step.gt.1d30) step=max(1d-3,time2-time1)

* Interpolate the orbit

      time=time1
100   continue
         i=getorb(time,lat,lon,orbit,path,.true.)
	 if (lon.gt.lon1) lon=lon-360d0
	 if (lon.lt.lon0) lon=lon+360d0
	 if (lat.ge.lat0 .and. lat.le.lat1 .and. lon.le.lon1) then
	    call geoxyz(lat*rad,lon*rad,orbit,xyz,r)
            write (*,600) time,i,lat,lon,orbit,xyz
	    skip=.false.
	 else if (.not.skip) then
	    write (*,601)
	    skip=.true.
         endif
	 time=time+step
      if (time.le.time2) goto 100
600   format (f14.3,i3,2f13.7,f12.3,3f13.3)
601   format ('>')

* End

9999  end
