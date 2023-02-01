*-----------------------------------------------------------------------
**I8SWAP -- Swap one or more 8-byte integers
*+
      SUBROUTINE I8SWAP (N, I)
      INTEGER*4 N
*     INTEGER*8 I(N)
*
* This routine swaps a 8-byte integer variable or the elements of a
* 8-byte integer array, i.e. it converts their representation between
* little and big endian (in either direction).
*
* Arguments:
*  N  (input) : Number of elements to swap
*  I  (input) : Integer or integer array to swap
*    (output) : Swapped integer or integer array
*
* Examples:
*     CALL I8SWAP(1,I)    ! Swap integer I
*     CALL I8SWAP(4,I)    ! Swap integer elements I(1) through I(4)
*     CALL I8SWAP(3,I(7)) ! Swap integer elements I(7) through I(9)
*-
*  9-Aug-2000 - Created by Remko Scharroo (DUT/DEOS)
*-----------------------------------------------------------------------
      INTEGER*1 I(8,N),B
      INTEGER*4 K
      DO K=1,N
         B=I(1,K)
         I(1,K)=I(8,K)
         I(8,K)=B
         B=I(2,K)
         I(2,K)=I(7,K)
         I(7,K)=B
         B=I(3,K)
         I(3,K)=I(6,K)
         I(6,K)=B
         B=I(4,K)
         I(4,K)=I(5,K)
         I(5,K)=B
      ENDDO
      END
