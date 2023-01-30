C    ******************************************************************
C                                                      44. MATMMM
C
C    PURPOSE
C    TO MULTIPLY TWO REAL MATRICES  AB=C
C
C    USAGE
C    CALL MATMMM(M,N,K,L,A,B,C)
C
C    DESCRIPTION OF PARAMETERS
C    M    I SECOND DIMENSION OF B AND C
C    N    I FIRST DIMENSION OF A AND C
C    K    I FIRST DIMENSION OF B AND SECOND DIMENSION OF A
C    L    I FIRST DIMENSION OF A, B AND C AS DECLARED IN THE CALLING
C           (SUB)PROGRAM,L AT LEAST M, N OR K
C    A    I REAL ARRAY OF DIMENSION (N,K)
C    B    I REAL ARRAY OF DIMENSION (K,M)
C    C    O RESULT, REAL ARRAY OF DIMENSION (N,M)
C
C    *******************************************************************
      SUBROUTINE MATMMM(M,N,K,L,A,B,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER R
      DIMENSION A(L,K),B(L,M),C(L,M)
      DO 80 I=1,N
         DO 80 J=1,M
         T=0D0
         DO 40 R=1,K
  40        T=T+A(I,R)*B(R,J)
  80     C(I,J)=T
      END
