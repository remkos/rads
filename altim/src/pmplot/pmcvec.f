**PMCVEC -- convert coordinates of several vectors.
*+
      SUBROUTINE PMCVEC (N, X, Y, A)
      INTEGER N
      REAL X(N), Y(N), A(N)
*
* Convert real-world coordinates (longitude,latitude) and the azimuth of a
* vector to map coordinates (x,y) and a tilt.
* The real-world coordinates are the position coordinates in degrees, and the
* azimuth is measured eastward from the north in degrees. The map coordinates
* are mapped linearly in the viewport and depend on the projection type used;
* the tilt is measured anti-clockwise (in degrees) starting horizontal pointing
* to the right. Even for linear projections (where the real-world coordinates
* are the same as the map coordinates) the azimuth is different from the tilt,
* not only because they are measured in different ways, but also simply because
* the linear rectangular projections do not preserve angles.
* This routine is very useful as a preparatory step to PGMARK.
*
* Arguments:
*  N    (input) : Number of points of which the coordinates must be converted.
*  X,Y  (input) : Longitude and latitude of the points (degrees).
*      (output) : Map coordinates.
*  A    (input) : Direction of the vector, measured from the north eastward
*                 (in degrees) = Azimuth.
*      (output) : Tilt (degrees) corresponding to the azimuth, measured
*                 anti-clockwise starting horizontal, pointing to the right.
*--
* 23-Jan-1991 - created [Remko Scharroo]
* 19-Mar-1992 - Standardize PMPLOT
*  7-Jul-1994 - All variables defined
* 11-Jul-1996 - Adjusted to PGPLOT 5.1
*-----------------------------------------------------------------------
      INCLUDE 'pgplot.inc'
      INCLUDE 'pmplot.inc'
      REAL LON(3), LAT(3)
      REAL SINA, COSA, EPS, FX, FY
      INTEGER I
*
      EPS=0.5
      FX=PGXLEN(PGID)/(PGXTRC(PGID)-PGXBLC(PGID))
      FY=PGYLEN(PGID)/(PGYTRC(PGID)-PGYBLC(PGID))
*
      DO I=1,N
         SINA=EPS*SIN(A(I)*RAD)
         COSA=EPS*COS(A(I)*RAD)*COS(Y(I)*RAD)
         LON(1)=X(I)-SINA
         LON(2)=X(I)
         LON(3)=X(I)+SINA
         LAT(1)=Y(I)-COSA
         LAT(2)=Y(I)
         LAT(3)=Y(I)+COSA
         CALL PMCONV(3,LON,LAT)
         X(I)=LON(2)
         Y(I)=LAT(2)
         A(I)=ATAN2((LAT(3)-LAT(1))*FY,(LON(3)-LON(1))*FX)/RAD
      ENDDO
      END
