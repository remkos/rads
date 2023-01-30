C    ******************************************************************
C                                                      45. MATMMV
C
C    PURPOSE
C    TO MULTIPLY A VECTOR BY A MATRIX, AX=B
C
C    USAGE
C    CALL MATMMV(M,N,L,A,X,B)
C
C    DESCRIPTION OF PARAMETERS
C    M    I NUMBER OF ELEMENTS IN B
C    N    I NUMBER OF ELEMENTS IN X
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING
C           (SUB)PROGRAM (L.GE.M)
C    A    I REAL ARRAY OF DIMENSION (L,K) WHERE K.GE.N
C    X    I GIVEN VECTOR, REAL VECTOR OF DIMENSION AT LEAST N
C    B    I RESULT, REAL VECTOR OF DIMENSION AT LEAST M
C
C    *******************************************************************
      SUBROUTINE MATMMV(M,N,L,A,X,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(L,N),X(N),B(M)
      DO 80 I=1,M
         T=0D0
         DO 40 J=1,N
  40        T=T+A(I,J)*X(J)
  80     B(I)=T
      END
