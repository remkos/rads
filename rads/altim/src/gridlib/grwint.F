#ifdef DOUBLE
#define Real REAL*8
#define GRWINT GRWINT8
#define BIG 1d35
#else
#define Real REAL*4
#define GRWINT GRWINT4
#define BIG 1e35
#endif
**GRWINT -- Interpolate a two-dimensional array using a weight function
*+
      FUNCTION GRWINT (GRID, X, Y, NX, NY, MX, SIGMA, CUTOFF)
      INTEGER*4 NX, NY, MX
      Real      GRWINT, GRID(MX,*), X, Y, SIGMA, CUTOFF
*
* This subroutine interpolates a two-dimensional regular grid using the
* exponential weight function
*     W = EXP ( -(D/SIGMA)**2 )
* where D is the distance from the point to be interpolated to any of
* the neighbouring grid points, SIGMA is a smoothing factor.
* The number of points that are considered 'neighbouring' is determined
* by the argument CUTOFF, which is the minimum weight of grid points to be
* used in the interpolation.
* The dimension of the grid is NX by NY. The real coordinates of the point to be
* interpolated are (X,Y) where 1 <= X <= NX and 1 <= Y <= NY.
*
* Arguments:
* GRID    (input): Two dimension array containing grid.
* X,Y     (input): Real grid coordinates of the point to be interpolated.
* NX,NY   (input): X- and Y- dimension of the grid.
* MX      (input): First dimension of the array as defined in the calling
*                  routine. MX >= NX.
* SIGMA   (input): Smoothing factor (as given above) in grid distances.
* CUTOFF  (input): Minimum weight of points to be included in the weighting.
* GRWINT (output): Interpolated value at point (X,Y).
*-
* 21-Jan-1991 - Created - Remko Scharroo
*-----------------------------------------------------------------------
      real*8  HORIZ,SIGMA2,SUM,SUMW,W,D2
      INTEGER KX,KX0,KX1,KY,KY0,KY1

      SIGMA2=SIGMA**2
      HORIZ=SQRT(-LOG(CUTOFF)*SIGMA2)

      KX0=MAX(1,MIN(INT(X-HORIZ),NX))
      KX1=MAX(1,MIN(INT(X+HORIZ),NX))
      KY0=MAX(1,MIN(INT(Y-HORIZ),NY))
      KY1=MAX(1,MIN(INT(Y+HORIZ),NY))

      SUM=0D0
      SUMW=0D0

      DO KX=KX0,KX1
         DO KY=KY0,KY1
            D2=(KX-X)**2+(KY-Y)**2
               W=EXP(-D2/SIGMA2)
            IF (W.GE.CUTOFF) THEN
               SUMW=SUMW+W
               SUM=SUM+W*GRID(KX,KY)
            ENDIF
         ENDDO
      ENDDO
      IF (SUMW.EQ.0D0) THEN
         GRWINT=BIG
      ELSE
         GRWINT=SUM/SUMW
      ENDIF
      END
