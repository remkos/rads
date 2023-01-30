**PMX -- determine map X coordinate
*+
      REAL FUNCTION PMX (Y, LON)
      REAL Y, LON
*
* Compute map x-coordinate of a point of given map y-coordinate and longitude.
*
* Arguments:
*  Y    (input) : Map y-coordinate.
*  LON  (input) : Real-world longitude (degrees).
*  PMX (output) : Map x-coordinate.
*--
*  9-Jan-1991 - created [Remko Scharroo]
* 14-Jan-1991 - implement conic and azimuthal projections.
* 30-Oct-1991 - implement more projections.
* 13-Jan-1992 - Standardize PMPLOT.
*  2-Apr-1993 - Include tilted rectangular projection.
* 30-Jun-1994 - Include polar projections.
*-----------------------------------------------------------------------
      INCLUDE 'pmplot.inc'
      real a,f

      if (ptype.ge.1 .and. ptype.lt.10) then
         PMX=lon
      else if (azmtal) then
         PMX=1e30
      else if (ptype.eq.21 .or. ptype.eq.22) then
         f=factl*(lon-loncen)
         PMX=-y*tan(f)
      else if (ptype.eq.31) then
         PMX=(lon-loncen)*cos(y*rad)
      else if (ptype.eq.32) then
         PMX=(lon-loncen)*sqrt(1-(y/90.)**2)
      else if (ptype.eq.33) then
         f=2*atan(y/factk)/rad
         PMX=(lon-loncen)*cos(factl*f)
      else if (ptype.eq.34) then
	 PMX=lon+y*factl
      else if (ptype.eq.41) then
	 if (y.gt.0 .neqv. abs(lon-loncen).gt.90) then
	    pmx=1e30
	 else
	    a=(lon-loncen)*rad
	    pmx=-y*tan(a)
	 endif
      else if (ptype.eq.42) then
	 if (y.lt.0 .neqv. abs(lon-loncen).gt.90) then
	    pmx=1e30
	 else
	    a=(lon-loncen)*rad
	    pmx=y*tan(a)
	 endif
      else
         pmx=1e30
      endif
      end
