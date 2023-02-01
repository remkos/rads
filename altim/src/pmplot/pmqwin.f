**PMQWIN -- query area of the selected plot
*+
      SUBROUTINE PMQWIN (LON0, LON1, LAT0, LAT1)
      REAL LON0, LON1, LAT0, LAT1
*
* Returns the (theoretical) limits of the area in the current plot.
*
* Arguments:
*  LON0 (output) : The most western longitude of the area (degrees).
*  LON1 (output) : The most eastern longitude of the area (degrees).
*  LAT0 (output) : The most southern latitude of the area (degrees).
*  LAT0 (output) : The most northern latitude of the area (degrees).
*
*--
* 28-Jun-1995 - Created.
*-----------------------------------------------------------------------
      INCLUDE 'pmplot.inc'

      if (PMOPEN.LT.4) then
	 CALL GRWARN('PMQWIN: Use PMWINDOW first')
	 RETURN
      endif
      lon0=lonmin
      lon1=lonmax
      lat0=latmin
      lat1=latmax
      end
