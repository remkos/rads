**MATSY1 - Fill an integer array of pointers for packed matrices
*+
      SUBROUTINE MATSY1 (N, H)
      INTEGER N, H(N)
*
* This subroutine fills an integer with pointers that can be used to
* find the index of elements of an UPPER packed symmetric matrix as
* used in LAPACK and the (old) NUMLIB.
*
* Example: A(.) is an upper packed symmetric matrix. The element
* (I,J) has index H(I)+J for I>J or H(J)+I for I<J
*
* Remark: H(I)=I*(I-1)/2
*-
* 29-Jul-2000 - Copied from Numlib
*-----------------------------------------------------------------------
      INTEGER I
      H(1)=0
      DO I=2,N
         H(I)=H(I-1)+I-1
      ENDDO
      END
