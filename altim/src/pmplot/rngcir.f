**RNGCIR -- compute coordinates of points at equal distance to location
*+
      SUBROUTINE RNGCIR (LON, LAT, RADIUS, N, AZI1, AZIN, X, Y)
      INTEGER N
      REAL LON, LAT, RADIUS, AZI1, AZIN, X(N), Y(N)
*
* Compute the coordinates of a number of points with a specified distance to a
* certain location. In other words, compute the coordinates of a range circle.
*
* Arguments:
*  LON    (input) : Longitude (degrees) of the center of the range circle.
*  LAT    (input) : Latitude (degrees) of that point.
*  RADIUS (input) : Radius of the first range circle (kilometers).
*  N      (input) : Number of points on the range circle. (N should be greater
*                   than 1)
*  AZI1   (input) : Lower azimuth boundary (degrees) measured from the north
*                   eastward.
*  AZIN   (input) : Upper azimuth boundary (degrees).
*  X     (output) : Array of longitude coordinates (degrees) of points
*                   1 through N.
*  Y     (output) : Array of latitude coordinates (degrees) of points
*                   1 through N.
*--
* 17-Jan-1991 - created [Remko Scharroo]
* 15-Mar-1991 - No more jumps from -180 to +180 degrees or v.v.
*  7-Jul-1994 - All variables defined
*-----------------------------------------------------------------------
      REAL SESF, SECF, X0, A, EPS, RAD, SINE, CESF, CECF
      INTEGER J

      IF (N.LE.1) RETURN
*
      RAD=ATAN(1.)/45
      EPS=RADIUS/6371.
      SINE=SIN(EPS)
      CECF=COS(EPS)*COS(LAT*RAD)
      CESF=COS(EPS)*SIN(LAT*RAD)
      SECF=SINE*COS(LAT*RAD)
      SESF=SINE*SIN(LAT*RAD)
      X0=0
      DO 100 J=1,N
         A=(AZI1+(J-1)*(AZIN-AZI1)/(N-1))*RAD
         X(J)=ATAN2(SINE*SIN(A),CECF-SESF*COS(A))/RAD+LON
         IF (J.GT.1 .AND. ABS(X(J)-X0).GT.180)
     .      X(J)=X(J)-SIGN(360.,X(J)-X0)
         X0=X(J)
  100    Y(J)=ASIN(CESF+SECF*COS(A))/RAD
      END
