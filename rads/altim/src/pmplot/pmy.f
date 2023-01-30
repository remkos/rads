**PMY -- determine map Y coordinate
*+
      REAL FUNCTION PMY (X, LAT)
      REAL X, LAT
*
* Compute map y-coordinate of a point of given map x-coordinate and latitude.
*
* Arguments:
*  X    (input) : Map x-coordinate.
*  LAT  (input) : Real-world latitude (degrees).
*  PMY (output) : Map y-coordinate.
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
      integer j

      if (ptype.eq.1 .or. ptype.eq.31) then
         PMY=lat
      else if (ptype.eq.2) then
         PMY=factk*sin(lat*rad)
      else if (ptype.eq.3 .or. ptype.eq.4) then
         PMY=factk*log(tan(qpi+lat*factl))
      else if (ptype.eq.5 .or. ptype.eq.33) then
         PMY=factk*tan(lat*rad/2)
      else if (ptype.eq.6) then
         PMY=factk*sin(tan(lat*rad)/factl)
      else if (azmtal) then
         PMY=1e30
      else if (ptype.eq.21 .or. ptype.eq.22) then
         a=factk-lat
         f=max(-1.,min(1.,x/a))
         f=asin(f)
         PMY=-a*cos(f)
      else if (ptype.eq.32) then
         a=lat*rad
         f=pi*sin(a)
         do j=1,10
            a=f-sin(a)
         enddo
         PMY=sin(a/2)*90.
      else if (ptype.eq.34) then
	 PMY=lat*factk
      else if (ptype.eq.41) then
	 f=(90.-lat)*factk
	 a=asin(x/f)
	 PMY=-f*cos(a)
      else if (ptype.eq.42) then
	 f=(lat+90.)*factk
	 a=asin(x/f)
	 PMY=f*cos(a)
      else
         PMY=1e30
      endif
      end
