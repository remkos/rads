**VNORM -- Norm of a vector
*+
      FUNCTION VNORM (V)
      REAL*8 VNORM, V(3)
*
* This function returns the norm (or length) of a 3-dimensional vector V.
*
* Arguments:
*  V      (input): A 3-dimensional vector
*  VNORM (output): The norm (or length) of V
*-
*  5-Jul-1993 - Created
*-----------------------------------------------------------------------
      vnorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      end
