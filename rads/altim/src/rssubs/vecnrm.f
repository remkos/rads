**VECNRM -- Normalize a 3D vector
*+
      SUBROUTINE VECNRM (X, V)
      REAL*8 X(3), V(3)
*
* This routine computes the normalized vector V of a
* 3-dimensional vector X. ( V = X / || X || )
*
* Arguments:
*  X    (input): 3-dimensional vector
*  V   (output): Normalized vector
*-
* 20-Dec-1996 - Created
*-----------------------------------------------------------------------
      real*8 n,vnorm
      n=vnorm(x)
      v(1)=x(1)/n
      v(2)=x(2)/n
      v(3)=x(3)/n
      end
