**GRENTIDE -- Give tidal height at given location according to Grenoble model
*+
      FUNCTION GRENTIDE (UTC, LAT, LON, NAME)
      REAL*8 GRENTIDE, UTC, LAT, LON
      CHARACTER*(*) NAME
*
* This function returns the total Ocean Tidal elevation for a number of
* given tidal constituents, using the Grenoble model.
* Long periodic equilibrium tides are added.
* You must set the environment variable WAVE_PATH to the name of the
* directory where the tide model can be found.
*
* Arguments:
*  UTC       (input): Epoch in UTC seconds since 1.0 Jan 1985.
*  LAT       (input): Geodetic latitude of the location in degrees.
*  LON       (input): East longitude of the location in degrees.
*  NAME      (input): Name of the model ("fes95.2", "load", or "fes95.2+load")
*  GRENTIDE (output): Total Ocean Tidal elevation in metres.
*-
* 11-Jul-1995: Generated from flather.f linking otide.f
*  4-Oct-1995: New call (includes NAME)
* 17-Feb-1997: Relaxing number of values needed for interpolation from 4 to 3.
* 10-Dec-1998: Relaxing number of values needed for interpolation from 3 to 1.
*-----------------------------------------------------------------------
      real*4	tide,slon,slat
      real*8	tlp,t,utc85
      parameter	(utc85=46066d0*86400d0)
      integer	istat
      logical   addlpeqmt
      save      addlpeqmt
      include	'common.h'

* Initialise on first call

      if (name.ne.model) then
	 model=name
	 addlpeqmt=(name.ne.'error' .and. name.ne.'load')
	 call init_data
      endif

* Convert longitude to the -180 to 180 range.

      slat=lat
      slon=lon

      if (slon.gt.180.) then
	 slon=slon-360.
      else if (slon.lt.-180.) then
	 slon=slon+360.
      endif

* Call Grenoble tide model

      t=(utc/86400d0+31046d0)/36525d0
      call grenoble(t,slon,slat,tide,istat)

* Return error value if istat=0 (land or out of area),
* or istat < 1 (not even on area boundaries) or |value| > 9m
* Otherwise, convert tidal elevation in centimetres to metres.

      if (istat.lt.1 .or. abs(tide).gt.900.) then
	 grentide=1d30
      else if (addlpeqmt) then
	 call lpeqmt(utc+utc85,lat,tlp)
	 grentide=(tide+tlp)/1d2
      else
	 grentide=tide/1d2
      endif

      end
