**ISNAN -- Checks if value is Not-A-Number
*+
      FUNCTION ISNAN (X)
      implicit none
      LOGICAL ISNAN
      REAL*8  X
*
* The function ISNAN returns TRUE is X is NAN (Not-A-Number)
*-
* 25-Jan-2004 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      ISNAN=X.NE.X
      END
