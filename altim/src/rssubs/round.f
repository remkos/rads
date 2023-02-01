**ROUND -- Round double float off to certain precision
*+
      FUNCTION ROUND (X, DELTA)
      REAL*8 ROUND, X, DELTA
*
* This function rounds the value X off to a precision DELTA.
* For example, ROUND (1.0153d0, 1d-2) produces 1.02d0
*
* This function only works if the parts on either side of the decimal
* point fit within an INTEGER*4
*
* Input arguments:
*  X    : Original value that needs rounding
*  DELTA: Precision of the returned value
*
* Return value:
*  ROUND: Value of X rounded to the nearest step DELTA
*-
*  3-Jun-2008 - Created by Remko Scharroo (Altimetrics LLC)
*-----------------------------------------------------------------------
      integer*4 i,j
      i = nint(x)
      j = nint((x-i)/delta)
      round = i + j * delta
      end
