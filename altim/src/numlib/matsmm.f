C    ******************************************************************
C                                                      46. MATSMM
C
C    PURPOSE
C    TO MULTIPLY A MATRIX BY A MATRIX  AB=C,
C    WHERE A IS A SYMMETRIC MATRIX IN VECTOR FORM
C
C    USAGE
C    CALL MATSMM(N,M,L,H,A,B,C)
C
C    DESCRIPTION OF PARAMETERS
C    N    I DIMENSION OF A
C    M    I SECOND DIMENSION OF B AND C
C    L    I FIRST DIMENSION OF C AS DECLARED IN THE CALLING
C           (SUB)PROGRAM (L.GE.N)
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH SUBROUTINE MATSY1
C    A    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX A. THE DIMENSION IS AT LEAST N*(N+1)/2
C           A CAN BE FORMED BY SUBROUTINE MATSYD
C    B    I GIVEN MATRIX OF DIMENSION (L,N)
C    C    O RESULTANT MATRIX OF DIMENSION (L,N)
C
C    *******************************************************************
      SUBROUTINE MATSMM(N,M,L,H,A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER H(N)
      REAL*8 A(*),B(L,M),C(L,M)
      DO 40 I=1,N
         DO 40 J=1,M
            T=0D0
            DO 10 K=1,I
  10           T=T+A(K+H(I))*B(K,J)
            DO 20 K=I+1,N
  20           T=T+A(I+H(K))*B(K,J)
  40        C(I,J)=T
      END
