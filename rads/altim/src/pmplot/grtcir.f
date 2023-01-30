**GRTCIR -- compute coordinates of points on a great circle
*+
      SUBROUTINE GRTCIR (LON1, LAT1, LON2, LAT2, N, X, Y)
      INTEGER N
      REAL LON1, LAT1, LON2, LAT2, X(N), Y(N)
*
* Compute the coordinates of a number of points on a great circle through
* two points.
*
* Arguments:
*  LON1   (input) : Longitude (degrees) of the first point.
*  LAT1   (input) : Latitude (degrees) of the first point.
*  LON2   (input) : Longitude (degrees) of the second point.
*  LAT2   (input) : Latitude (degrees) of the second point.
*  N      (input) : Number of points on the great circle. (N should be greater
*                   than 1)
*  X     (output) : Array of longitude coordinates (degrees) of points
*                   1 through N.
*  Y     (output) : Array of latitude coordinates (degrees) of points
*                   1 through N.
*--
* 17-Jan-1991 - created [Remko Scharroo]
* 15-Mar-1991 - No more jumps from -180 to +180 degrees or v.v.
*-----------------------------------------------------------------------
      REAL X1, X2, Y1, Y2, Z1, Z2, RAD, COSE, EPS, SINE, C1, C2, C3, A
      INTEGER J

      IF (N.LE.1) RETURN
*
      RAD=ATAN(1.)/45
      X1=COS(LAT1*RAD)*COS(LON1*RAD)
      Y1=COS(LAT1*RAD)*SIN(LON1*RAD)
      Z1=SIN(LAT1*RAD)
      X2=COS(LAT2*RAD)*COS(LON2*RAD)
      Y2=COS(LAT2*RAD)*SIN(LON2*RAD)
      Z2=SIN(LAT2*RAD)
      COSE=X1*X2+Y1*Y2+Z1*Z2
      EPS=ACOS(COSE)
      SINE=SIN(EPS)
      C1=(X2-X1*COSE)/SINE
      C2=(Y2-Y1*COSE)/SINE
      C3=(Z2-Z1*COSE)/SINE
      DO J=1,N
         A=(J-1)*EPS/(N-1)
         X(J)=ATAN2(C2*SIN(A)+Y1*COS(A),C1*SIN(A)+X1*COS(A))/RAD
         IF (J.GT.1 .AND. ABS(X(J)-X(1)).GT.180)
     .       X(J)=X(J)-SIGN(360.,X(J)-X(1))
         Y(J)=ASIN(C3*SIN(A)+Z1*COS(A))/RAD
      ENDDO
      END
