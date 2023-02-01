**IRITEC -- Compute TEC from IRI model at 800 and 1400 km
*+
      SUBROUTINE IRITEC (UTC, LAT, LON, TEC1, TEC2, VERBOSE, TYPE)
      implicit none
      REAL*8	UTC, LAT, LON, TEC1, TEC2
      LOGICAL*4 VERBOSE
      CHARACTER*(*) TYPE

* This routine computes the total electron content based on an
* IRI model integrated up to 800 and 1400 km altitude.
* Input parameters are the geographical latitude and longitude
* and the time in seconds since 1 Jan 1985.
*
* Input arguments:
*   UTC    : UTC time in seconds since 1 Jan 1985
*   LAT    : Geographical latitude (deg)
*   LON    : Geographical longitude (deg)
*   VERBOSE: If .true. process info is provided
*   TYPE   : IRI model used (has to be "iri95")
*
* Output arguments:
*   TEC1 : Total electron content until  800 km (TEC units)
*   TEC2 : Total electron content until 1400 km (TEC units)
*
* (1 TEC unit = 10^16/m^2)
*-
* 18-Sep-2006 - Using NetCDF files
* 26-May-2004 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4	mjd0,mjd1,mjd2,mjd3,mjd,yy,mm,dd,map,kx,ky,i,j,iriread
      integer*2	tecmap(73,73,2,12,2)
      real*8	xmjd,dmonth,xmap,dmap,x,y,wul,wur,wll,wlr,
     |		nan,tec(2,2,2),units
      logical*4 exist
      save	mjd0,mjd1,nan,exist,units,tecmap
      data	mjd0/99999/,mjd1/-99999/,nan/-1d30/,exist/.false./

* Initialize whenever UTC is out of current range

      xmjd=utc/86400d0+46066d0
      mjd=int(xmjd)
      if (mjd.lt.mjd0 .or. mjd.gt.mjd1) then

* Compute mid-month boundaries (mjd2 and mjd3) for current epoch

	 call mjd2ymd(mjd,yy,mm,dd)
	 call ymd2mjd(yy,mm,15,mjd2)
	 if (mm.eq.2) mjd2=mjd2-1
	 if (mjd.lt.mjd2) then
	    mjd3=mjd2
	    call ymd2mjd(yy,mm-1,15,mjd2)
	    if (mm-1.eq.2) mjd2=mjd2-1
	 else
	    call ymd2mjd(yy,mm+1,15,mjd3)
	    if (mm+1.eq.2) mjd3=mjd3-1
	 endif

* Read the IRI95 maps for the new boundary epochs

	 exist=(iriread(mjd2,tecmap(1,1,1,1,1),units,verbose,type).eq.0
     |	  .and. iriread(mjd3,tecmap(1,1,1,1,2),units,verbose,type).eq.0)
	 mjd0=mjd2
	 mjd1=mjd3
	 nan=sqrt(nan)
      endif

* If either monthly grid does not exist, return NaN

      if (.not.exist) then
         tec1=nan
	 tec2=nan
	 return
      endif

* Compute the index for the first hourly map for interpolation

      xmap=(xmjd-mjd)*12+1
      map=int(xmap)
      dmap=xmap-map
      dmonth=(xmjd-mjd0)/(mjd1-mjd0)

* Compute grid indices within the first hourly map (after rotation)

      x=lon+dmap*30d0
      if (x.lt.-180d0) x=x+360d0
      if (x.ge.+180d0) x=x-360d0
      x=(x+180d0)/5d0+1
      kx=int(x)
      y=(lat+90d0)/2.5d0+1
      ky=int(y)
      x=x-kx
      y=y-ky
      wll=(1-x)*(1-y)
      wlr=   x *(1-y)
      wul=(1-x)*   y
      wur=   x *   y
*      write (*,*) kx,x,ky,y,map,dmap,dmonth

* Interpolate the TEC in the first hourly map (both levels and both months)

      do j=1,2
	 do i=1,2
	    tec(1,i,j)=
     |		wll*tecmap(kx,ky  ,i,map,j)+wlr*tecmap(kx+1,ky  ,i,map,j)+
     |		wul*tecmap(kx,ky+1,i,map,j)+wur*tecmap(kx+1,ky+1,i,map,j)
	 enddo
      enddo

* Rotate the X coordinate in the second hourly map (both levels and both months)

      kx=kx-6
      if (kx.lt.1) kx=kx+72
      map=map+1
      if (map.gt.12) map=map-12
      do j=1,2
	 do i=1,2
	    tec(2,i,j)=
     |		wll*tecmap(kx,ky  ,i,map,j)+wlr*tecmap(kx+1,ky  ,i,map,j)+
     |		wul*tecmap(kx,ky+1,i,map,j)+wur*tecmap(kx+1,ky+1,i,map,j)
	 enddo
      enddo

* Interpolate in time (between hourly maps and months) for both levels

      tec1=((tec(1,1,1)*(1-dmap)+tec(2,1,1)*dmap)*(1-dmonth)+
     |      (tec(1,1,2)*(1-dmap)+tec(2,1,2)*dmap)*dmonth)*units
      tec2=((tec(1,2,1)*(1-dmap)+tec(2,2,1)*dmap)*(1-dmonth)+
     |      (tec(1,2,2)*(1-dmap)+tec(2,2,2)*dmap)*dmonth)*units
      end

**IRIREAD -- Read IRI Ionosphere Map
*+
      FUNCTION IRIREAD (MJD, TECMAP, UNITS, VERBOSE, TYPE)
      INTEGER*4 MJD, IRIREAD
      INTEGER*2 TECMAP(*)
      REAl*8	UNITS
      LOGICAL	VERBOSE
      CHARACTER TYPE*(*)

* This routine loads a TEC map stored in an IRI file into memory.
* If the requested file does not exist or can not be read, IRIREAD
* returns a non-zero value, and the first element in TECMAP is
* set to -9999.
*
* It is assumed that all grids cover the same area and
* have the same grid spacing and the each daily file contains
* twelve TEC maps.
*
* Longitude: -180 to +180 with intervals of 5 degrees (73 grid lines)
* Latitude:  -90.0 to +90.0 with intervals of 2.5 degrees (73 grid lines)
* TEC:       In units of "UNITS" TEC Units, or  UNITS*10**16 electrons/m**2
*
* IRI files should be stored in $ALTIM/data/TYPE/maps, where ALTIM
* is an environment variable containing a directory name chosen by
* the user, with /user/altim as default. TYPE currently has to be
* "iri95", but other types may be considered later.
* The IRI files are named: iri_YYMMDD.nc, where YY is the year,
* MM the month and DD the day number.
*
* Arguments:
*  MJD      (input) : Modified Julian Date of IRI file
*  TECMAP  (output) : Integer array with the 12 TEC maps
*  UNITS   (output) : Units of the TECMAP array (in TEC units, usually 0.1)
*  VERBOSE  (input) : If .true., be verbose
*  TYPE     (input) : Source of IRI file (has to be "iri95")
*  IRIREAD (output) : Status report: 0 = no error, 1 = error opening file,
*                     2 = error reading file.
*-----------------------------------------------------------------------
      integer*4 lnblnk,l,i,yy,mm,dd,ncid,varid
      character*256 filenm
      include "netcdf.inc"

* Construct the pathname of the IRI file

      filenm='/user/altim'
      call checkenv('ALTIM',filenm,l)
      i=lnblnk(type)
      call mjd2ymd(mjd,yy,mm,dd)
      write (filenm(l+1:),600) type(:i),mod(yy,100),mm,dd
600   format ('/data/',a,'/maps/iri_',3i2.2,'.nc')
610   format ('(IRITEC: Reading ',a,' ... ',$)
620   format (a,')')
      if (verbose) write (0,610) filenm(:l+i+25)

* Open the file and load the grid

      if (nf_open(filenm,nf_nowrite,ncid).ne.nf_noerr) then
         if (verbose) write (0,620) 'error opening file'
	 tecmap(1)=-9999
	 iriread=1
	 return
      endif
      if (nf_inq_varid(ncid,"tec",varid).ne.nf_noerr .or.
     |	nf_get_var_int2(ncid,varid,tecmap).ne.nf_noerr .or.
     |	nf_get_att_double(ncid,varid,"scale_factor",units).ne.nf_noerr)
     |	then
         if (verbose) write (0,620) 'error reading file'
	 tecmap(1)=-9999
	 iriread=2
	 return
      endif
      i=nf_close(ncid)

* Swap bytes if necessary

      iriread=0
      if (verbose) write (0,620) 'done'
      end
