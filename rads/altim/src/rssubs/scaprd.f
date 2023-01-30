**SCAPRD -- Scalar product of two 3D vectors
*+
      FUNCTION SCAPRD (X, Y)
      REAL*8 SCAPRD, X(3), Y(3)
*
* This function returns the scalar (or inner) product of two 3-dimensional
* vectors X and Y
*
* Arguments:
*  X, Y    (input): Two 3-dimensional vectors
*  SCAPRD (output): Scalar (or inner) product of X and Y (X.Y)
*-
*  5-Jul-1993 - Created
*-----------------------------------------------------------------------
      scaprd=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
      end
