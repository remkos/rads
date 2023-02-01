**GIMTEC -- Compute TEC from Global Ionosphere Model
*+
      FUNCTION GIMTEC (UTC, LAT, LON, NDAYS, TECMAP, VERBOSE, TYPE)
      INTEGER*4 NDAYS
      REAL*8    UTC, LAT, LON, GIMTEC
      INTEGER*2 TECMAP(73,71,NDAYS*12+1)
      LOGICAL   VERBOSE
      CHARACTER TYPE*4

* Compute the Total Electron Content (TEC) based on the GPS Ionosphere
* Model maps of TEC. This routine loads the appropriate IONEX files from
* $ALTIM/data/gim and makes a spatio-temporal interpolation in the
* two closest TEC maps, using the solar-geographical longitude as
* longitude coordinate. The TEC maps are spaced by 2 hours.
*
* Three sources of IONEX files are supported, specified by the argument
* TYPE: "JPLG" for JPL analysis grids, "JPLQ" for JPL quick-look grids,
* "CODG" for CODE analysis grids.
*
* The routine will search for the IONEX files in the directory
* $ALTIM/data/gim, where $ALTIM is an environment variable. If this
* variable is not set, the default /user/altim is used. Files are stored
* is yearly directories, $ALTIM/data/gim/YYYY, where YYYY is the year
* (e.g. 2002). File names are:
* - JPLGddd0.yyI.gz (JPL analysis files, Gnu compressed)
* - JPLQddd0.yyI.gz (JPL quick-look files, Gnu compressed)
* - CODGddd0.yyI.Z  (CODE analysis files, Unix compressed)
* where yy is the year (e.g. 02) and ddd the day number (e.g. 011)
*
* Input to the routine are the satellite geographical latitude (LAT),
* longitude (LON) and the epoch of observation in seconds since 1 Jan
* 1985 (UTC). Also the data source TYPE has to be specified (see above).
*
* Returned is the TEC between the satellite and its nadir. When either
* of the two TEC maps is not available for interpolation, a value of
* 1D30 is returned. To convert this TEC into an ionospheric path delay
* (positive, in metres), compute:
*     IONO = C * TEC / F**2
* where F is the frequency of the signal in Hz and C is 40250d13
*
* To improve efficiency, more than one day can be loaded into memory.
* Choose NDAYS such that GIM files have to be loaded infrequently.
* To store the GIMs provide an array TECMAP with at least
* 73 * 71 * (NDAYS*12 + 1) 2-byte integers.
*
* This routine is loosely based on routines created by AIUB. See:
* http://www.cx.unibe.ch/aiub/ionosphere.html
*
* This routine runs efficiently only when time progresses forward.
* Temporary uncompressed files are written to /tmp.
* Execution requires the "gunzip" program.
* Files are opened on a self-chosen free unit number.
*
* Input arguments:
*   UTC    : Time in seconds since 1 Jan 1985
*   LAT    : Geographical latitude in degrees
*   LON    : Geographical longitude in degrees
*   NDAYS  : Number of days to be loaded into memory
*   TECMAP : Memory buffer to contain all TEC maps
*   VERBOSE: If .true. process info is provided
*   TYPE   : Source of IONEX files: "JPLG", "JPLQ", or "CODG"
*
* Output arguments:
*   GIMTEC : Total Electron Content between the satellite and its
*            nadir in TEC Units (10**16 electrons/m**2).
*-
* 23-May-2002 - Created by Remko Scharroo
* 16-Nov-2002 - Allow for reading of 13 maps per file
* 19-Nov-2002 - Small tweak to avoid discontinuities
*  4-Dec-2002 - Added JPLQ. Updated manual.
*  4-Feb-2003 - Corrected bug in dimension of oldtype
* 12-Mar-2003 - Allow use of special file for 2 Nov 2002
* 27-May-2004 - Set undetermined TEC to NaN
* 24-Jan-2008 - Introduced "diagonal" interpolation to better follow
*               geomagnetic lines
*-----------------------------------------------------------------------
      integer*4 mjd,map,kx,ky,gimread,mdate,k,hour1,mjd13
      real*8    x,y,xmap,dmap,dx,dy,wll,wlr,wul,wur,tec1,tec2,utc13,upy
      real*8    utc1/1d30/,utc2/-1d30/,dtec,nan/-1d30/,rad,lonp
      character oldtype*4/'----'/
      parameter (dtec=1d-1,mjd13=52580-46066,utc13=mjd13*86400d0,
     |		rad=1.74532925d-2,lonp=288.63d0)
      save      utc1,utc2,oldtype,nan

* Check if the maps corresponding to the time are loaded.
* - utc1 = UTC seconds of the first map loaded
* - utc2 = UTC seconds of the last map loaded

      if (utc.le.utc1 .or. utc.ge.utc2 .or. type.ne.oldtype) then
      
* Load as many grids as possible in the provided buffer space.
* Start with the day of the measurement.
* Make provisions for both 12 and 13 maps per day, and for the transition
* on MJD 52580 (2 Nov 2002)

	 if (utc.ge.utc13) then
	    hour1=0
	 else
	    hour1=1
	 endif
	 mjd=int((utc-hour1*3600d0)/86400d0)
	 utc1=mjd*86400d0+hour1*3600d0
	 if (verbose) write (0,600) type
	 do map=1,ndays*12,12
	    if (verbose) write (0,601) mdate(1,mjd+46066)
	    if (gimread(mjd+46066,type,tecmap(1,1,map)).ne.0) then
	       do k=0,12
	          tecmap(1,1,map+k)=-9999
	       enddo
	    endif
	    mjd=mjd+1
	    if (mjd.eq.mjd13) then
	       hour1=2
	       goto 100
	    endif
	 enddo
100	 oldtype=type
	 utc2=mjd*86400d0+(hour1-2)*3600d0
	 nan=sqrt(nan)
	 if (verbose) write (0,602)
      endif

* Compute the index for the first map for interpolation

      xmap=(utc-utc1)/7200d0+1
      map=int(xmap)
      dmap=xmap-map
      if (tecmap(1,1,map).eq.-9999 .or. tecmap(1,1,map+1).eq.-9999) then
         gimtec=nan
	 return
      endif

* Compute grid indices within the first map (after rotation)

      upy = 30d0 * 0.2d0 * sin((lon-lonp)*rad) * cos(lat*rad)
      x=lon+dmap*30d0
      if (x.lt.-180d0) x=x+360d0
      if (x.gt.+180d0) x=x-360d0
      x=(x+180d0)/5d0+1
      kx=int(x)
      dx=x-kx
      y=(lat+upy*dmap+87.5d0)/2.5d0+1
      ky=int(y)
      dy=y-ky
      wll=(1-dx)*(1-dy)
      wlr=   dx *(1-dy)
      wul=(1-dx)*   dy
      wur=   dx *   dy

* Interpolate the TEC in the first map

      tec1=wll*tecmap(kx,ky  ,map)+wlr*tecmap(kx+1,ky  ,map)+
     |     wul*tecmap(kx,ky+1,map)+wur*tecmap(kx+1,ky+1,map)
      tec1=tec1*dtec

* Rotate the X and Y coordinate in the second map and interpolate

      kx=kx-6
      if (kx.lt.1) kx=kx+72
      y=(lat-upy*(1d0-dmap)+87.5d0)/2.5d0+1
      ky=int(y)
      dy=y-ky
      wll=(1-dx)*(1-dy)
      wlr=   dx *(1-dy)
      wul=(1-dx)*   dy
      wur=   dx *   dy

      map=map+1
      tec2=wll*tecmap(kx,ky  ,map)+wlr*tecmap(kx+1,ky  ,map)+
     |     wul*tecmap(kx,ky+1,map)+wur*tecmap(kx+1,ky+1,map)
      tec2=tec2*dtec

* Interpolate in time between the two values

      gimtec=tec1*(1-dmap)+tec2*dmap

* Formats

600   format ('(Loading ',a,' TEC maps for day',$)
601   format (1x,i6.6,$)
602   format (')')
      end
