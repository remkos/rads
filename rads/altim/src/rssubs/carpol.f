**CARPOL -- Convert Cartesian to Polar
*+
      SUBROUTINE CARPOL (VECTOR, LAT, LON, R)
      REAL*8 VECTOR(3), LAT, LON, R
*
* Convert Cartesian VECTOR() = (X,Y,Z) to polar coordinates (LAT,LON,R)
*
* Arguments:
* VECTOR  (input) : Vector of X, Y, and Z coordinates.
* LAT    (output) : Latitude in radians
*                   = angle from XY-plane towards +Z-axis
* LON    (output) : Longitude in radians
*                   = angle in XY-plane measured from +X-axis towards +Y-axis
* R      (output) : Radius
*-
* 28-Apr-1991 - Created: Remko Scharroo, DUT/SOM (c)
*-----------------------------------------------------------------------
      R=SQRT(VECTOR(1)**2+VECTOR(2)**2+VECTOR(3)**2)
      LON=ATAN2(VECTOR(2),VECTOR(1))
      LAT=ASIN(VECTOR(3)/R)
      END
