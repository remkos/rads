**THETOP -- Computes the top of a parabola
*+
      SUBROUTINE THETOP (X0, X1, X2, T, XT)
      REAL*8 X0, X1, X2, T, XT
*
* Given three equidistant points, this routine computes the top of a
* parabola fitting through these points.
* The three input points are:
*     X0 = X(t=-1) ; X1 = X(t= 0) ; X2 = X(t=+1)
* The output will be the time T and the value XT of the top of the parabola,
* where time refers to the mid point X1 and is measured in time intervals
* between X0 and X1 (or X1 and X2). So  XT = X(t=T)
*
* This routine is very handy to determine if points are an ascending or
* descending track. E.G., if you have latitudes LAT(1...N), then point I
* is on an ascending track if ASCEN=.TRUE., after
*
*     CALL THETOP (LAT(I-1), LAT(I), LAT(I+1), T, XT)
*     ASCEN = (T*(LAT(I)-XT).LT.0)
*
* Arguments:
*   X0  (input): first input point  = X(t=-1)
*   X1  (input): second input point = X(t= 0)
*   X2  (input): third input point  = X(t=+1)
*   T  (output): time of the extreme
*   XT (output): value at the extreme = X(t=T)
*-
* 23-Jul-93 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      real*8 b,c
      b=(x2-x0)/2
      c=x2-b-x1
      t=-b/(2*c)
      xt=x1+b*t+c*t*t
      end
