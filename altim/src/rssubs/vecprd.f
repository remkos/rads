**VECPRD -- Vectorial product of two 3D vectors
*+
      SUBROUTINE VECPRD (X, Y, V)
      REAL*8 X(3), Y(3), V(3)
*
* This routine computes the vectorial (or outer) product V of two
* 3-dimensional vectors X and Y (V = X * Y)
*
* Arguments:
*  X, Y (input): Two 3-dimensional vectors
*  V   (output): Vectorial (or outer) product of X and Y
*-
*  5-Jul-1993 - Created
*-----------------------------------------------------------------------
      v(1)=x(2)*y(3)-x(3)*y(2)
      v(2)=x(3)*y(1)-x(1)*y(3)
      v(3)=x(1)*y(2)-x(2)*y(1)
      end
