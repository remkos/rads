**I4SWAP -- Swap one or more 4-byte integers
*+
      SUBROUTINE I4SWAP (N, I)
      INTEGER*4 N
*     INTEGER*4 I(N)
*
* This routine swaps a 4-byte integer variable or the elements of a
* 4-byte integer array, i.e. it converts their representation between
* little and big endian (in either direction).
*
* Arguments:
*  N  (input) : Number of elements to swap
*  I  (input) : Integer or integer array to swap
*    (output) : Swapped integer or integer array
*
* Examples:
*     CALL I4SWAP(1,I)    ! Swap integer I
*     CALL I4SWAP(4,I)    ! Swap integer elements I(1) through I(4)
*     CALL I4SWAP(3,I(7)) ! Swap integer elements I(7) through I(9)
*-
*    Aug-1991 - Created: Remko Scharroo, DUT/SSR&T (c)
*  2-Jul-1992 - Documented version
* 12-Nov-1993 - CHARACTER*1 replaced by BYTE
*-----------------------------------------------------------------------
      INTEGER*1 I(4,N),B
      INTEGER*4 K
      DO K=1,N
         B=I(1,K)
         I(1,K)=I(4,K)
         I(4,K)=B
         B=I(2,K)
         I(2,K)=I(3,K)
         I(3,K)=B
      ENDDO
      END
