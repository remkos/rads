**PGNORM - normalize X and Y coordinates
*+
      SUBROUTINE PGNORM (R, X, Y)
      REAL R, X, Y
*
* This routine normalizes the X and Y coordinates of a vector (X,Y),
* such that it obtains a length R.
*
* Arguments:
*  R           (input): New length of the vector.
*  X, Y (input/output): Coordinates of the vector.
*-
* 18-Mar-1992 - Created [Remko Scharroo].
*-----------------------------------------------------------------------
      REAL FACT
      FACT=R/SQRT(X*X+Y*Y)
      X=FACT*X
      Y=FACT*Y
      END
