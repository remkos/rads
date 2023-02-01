**GROUNDTRACK -- Compute the groundtrack positions of a satellite
*+
      program groundtrack
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
      implicit none

      real*8 pi,rad
      real*8 start,end,step,a,e,incl,og0,ogdot,w0,wdot,wm0,wmdot
      real*8 og,w,wm,theta,anmtot,epoch
      real*8 p,t,lat,lon,r,h,v(3)

      integer nstep,n

      pi=4*atan(1d0)
      rad=pi/180

      read (5,*) start,end,step
      read (5,*) epoch
      read (5,*) a
      read (5,*) e
      read (5,*) incl
      read (5,*) og0
      read (5,*) ogdot
      read (5,*) w0
      read (5,*) wdot
      read (5,*) wm0
      read (5,*) wmdot

      p=a*(1-e*e)
      w=w*rad
      incl=incl*rad
      step=step/86400
      nstep=nint((end-start)/step)

      do n=0,nstep
	 t=start+n*step
	 og=(og0+ogdot*(t-epoch))*rad
	 w=(w0+wdot*(t-epoch))*rad
	 wm=(wm0+wmdot*(t-epoch))*rad
	 theta=anmtot(wm-w,e)
         v(1)=p/(1+e*cos(theta))
	 v(2)=0
	 v(3)=0
	 call rotate(3,-(w+theta),v,v)
	 call rotate(1,-incl,v,v)
	 call rotate(3,-og,v,v)
	 call xyzgeo(v,r,lat,lon,h)
	 lat=lat/rad
	 lon=lon/rad
	 if (lon.lt.0) lon=lon+360
	 write (6,600) t,lat,lon,h
      enddo
  600 format (3f12.6,f12.3)
      end
