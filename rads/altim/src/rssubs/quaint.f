**QUAINT -- Quadratic interpolation
*+
      FUNCTION QUAINT (T, H, TT)
      REAL*8 QUAINT, T(0:2), H(0:2), TT
*
* Produce interpolated value QUAINT on epoch TT from tables H() and T(),
* using parabolic fit through 3 points (need not be equidistant).
*
* Arguments:
*  T      (input): Time tags for values in table H
*  H      (input): Values on time tags T
*  TT     (input): time-tag on which interpolated value must be determined
*  QUAINT (output): Interpolated value.
*-
*  5-Jul-1993 - New manual.
*-----------------------------------------------------------------------
      real*8 a,b,c
      a=(tt-t(0))/(t(2)-t(1))
      b=(tt-t(1))/(t(2)-t(0))*(h(2)-h(0))
      c=(tt-t(2))/(t(1)-t(0))*(h(1)-h(0))
      quaint=h(0)+a*(b-c)
      end
