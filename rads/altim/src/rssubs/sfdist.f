**SFDIST -- Spherical distance between two points
*+
      FUNCTION SFDIST (LAT0, LON0, LAT1, LON1)
      REAL*8 LAT0, LON0, LAT1, LON1, SFDIST
*
* Compute spherical distance between two points of given geocentric
* latitude and longitude.
* To avoid rounding errors at small distances, we use the haversine
* formula which is optimized for computing very small spherical
* distances, but only has inaccuracies determining distances to
* antipodal points.
*
* References:
* - R.W. Sinnott, "Virtues of the Haversine", Sky and Telescope,
*   vol. 68, no. 2, p. 159, 1984.
* - http://en.wikipedia.org/wiki/Great-circle_distance
*
* Arguments:
*  LAT0, LON0  (input): Geocentic latitude and longitude of one point (rad).
*  LAT1, LON1 (output): Geocentic latitude and longitude of other point (rad).
*  SFDIST     (output): Spherical distance in radians.
*-
* 31-Oct-1990: Created [Remko Scharroo].
* 15-Jun-2012: Bug removed in implementation of haversine formula.
*-----------------------------------------------------------------------
      REAL*8 S
      S=DSIN((LAT0-LAT1)/2)**2
     | +DCOS(LAT0)*DCOS(LAT1)*DSIN((LON0-LON1)/2)**2
      SFDIST=2*DASIN(DSQRT(S))
      END
