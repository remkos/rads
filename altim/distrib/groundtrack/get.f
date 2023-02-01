      real*8 lat,lon,time,orbit
      integer i,getorb
      logical verbose
      character*80 path

1     read (5,*) time
      time=(time-46066)*86400
      call getorb2(time,lat,lon,orbit,'ers1.elem',.true.)
      path='/home/ers1/ODR/JGM-2'
      i=getorb(time,lat,lon,orbit,path,verbose)
      write (6,600) time,lat,lon,orbit
  600 format (3f12.6,f12.3)
      goto 1
      end

**GETORB2 -- Simulate the groundtrack positions of a satellite
*+
      SUBROUTINE GETORB2 (TIME, LAT, LON, ORBIT, PATH, VERBOSE)
      REAL*8 TIME, LAT, LON, ORBIT
      CHARACTER*(*) PATH
      LOGICAL VERBOSE
*
* This program determines the goundtrack of a satellite, given a number
* of orbital parameters. These parameters are to be read from standard
* input. An input file could look like this:
*
* 49146.000000 49151.500000   300.0 <- period (MJD), step (sec)
*   49146.000000 <- epoch (MJD)
* 7163807.993057 <- semi-major axis (m)
*  .001161259861 <- eccentricity (-)
*   98.558372617 <- inclination (deg)
*  -22.499587526 <- node at epoch (deg), wrt Greenwich
* -359.999242458 <- node rate (deg/day), wrt Greenwich
*   91.000632829 <- perigee at epoch (deg)
*     .000000000 <- perigee rate (deg/day)
*   87.535233012 <- mean longitude at epoch (deg)
* 5153.123608154 <- angular motion (deg/day)
*
* The output of the program is produced at standard output and contains:
*
* Time            Latitude   Longitude    Altitude
* 49146.000000   81.146266  263.682797  798241.518
* 49146.003472   72.477277  184.532009  797047.514
* 49146.006944   55.811383  167.722291  793240.565
* ....
* 49151.496528  -28.092815  163.334700  794447.976
* 49151.500000  -10.358527  159.071475  787994.586
*
* Time is in MJD, the positions refer to the GRS80 ellipsoid
* (ae = 6378.137 km, 1/f = 298.257)
*-
*  6-Jul-1993 - Created by Remko Scharroo (DUT/SSR&T)
*-----------------------------------------------------------------------
      include "math.inc"

      integer n,i,narc,jarc
      parameter (n=999)
      real*8 start(n),end(n),a(n),e(n),incl(n),og0(n),ogdot(n),
     .		w0(n),wdot(n),wm0(n),wmdot(n),epoch(n)
      real*8 og,w,wm,theta,anmtot
      real*8 r,v(3)
      character*160 oldpath/' '/

      integer mjd85
      parameter (mjd85=46066)

      if (oldpath.ne.path) then
         open(35,file=path,status='old')
	 oldpath=path
	 narc=0
10       narc=narc+1
         read (35,*,end=90) start(narc),end(narc)
         read (35,*) epoch(narc)
         read (35,*) a(narc)
         read (35,*) e(narc)
         read (35,*) incl(narc)
         read (35,*) og0(narc)
         read (35,*) ogdot(narc)
         read (35,*) w0(narc)
         read (35,*) wdot(narc)
         read (35,*) wm0(narc)
         read (35,*) wmdot(narc)
	 read (35,*)
	 read (35,*)
         start(narc)=(start(narc)-mjd85)*86400
         end(narc)=(end(narc)-mjd85)*86400
         epoch(narc)=(epoch(narc)-mjd85)*86400
         a(narc)=a(narc)*(1-e(narc)*e(narc))
         incl(narc)=incl(narc)*rad
         og0(narc)=og0(narc)*rad
         ogdot(narc)=ogdot(narc)*rad/86400
         w0(narc)=w0(narc)*rad
         wdot(narc)=wdot(narc)*rad/86400
         wm0(narc)=wm0(narc)*rad
         wmdot(narc)=wmdot(narc)*rad/86400
	 goto 10
90       narc=narc-1
         close (35)
         jarc=1
      endif

100   if (time.lt.start(jarc)) then
	 do i=jarc,1,-1
	    if (time.ge.start(i) .and. time.le.end(i)) then
	       jarc=i
	       goto 110
	    endif
	 enddo
	 stop 'getorb2: no orbital elements for that epoch'
	 if (verbose) write (0,*) 'using arc ',jarc
      else if (time.gt.end(jarc)) then
	 do i=jarc,narc,1
	    if (time.ge.start(i) .and. time.le.end(i)) then
	       jarc=i
	       goto 110
	    endif
	 enddo
	 stop 'getorb2: no orbital elements for that epoch'
	 if (verbose) write (0,*) 'using arc ',jarc
      endif

110   continue

      og=og0(jarc)+ogdot(jarc)*(time-epoch(jarc))
      w=w0(jarc)+wdot(jarc)*(time-epoch(jarc))
      wm=wm0(jarc)+wmdot(jarc)*(time-epoch(jarc))
      theta=anmtot(wm-w,e(jarc))
      v(1)=a(jarc)/(1+e(jarc)*cos(theta))
      v(2)=0
      v(3)=0
      call rotate(3,-(w+theta),v,v)
      call rotate(1,-incl(jarc),v,v)
      call rotate(3,-og,v,v)
      call xyzgeo(v,r,lat,lon,orbit)
      lat=lat/rad
      lon=lon/rad
      write (6,600) time,lat,lon,orbit
  600 format (3f12.6,f12.3)
      end
