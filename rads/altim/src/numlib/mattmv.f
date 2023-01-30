C    ******************************************************************
C                                                     45B. MATTMV
C
C    PURPOSE                                       T
C    TO MULTIPLY A VECTOR BY A TRANSPOSED MATRIX, A X=B
C
C    USAGE
C    CALL MATTMV(M,N,L,A,X,B)
C
C    DESCRIPTION OF PARAMETERS
C    M    I NUMBER OF ELEMENTS IN B = SECOND DIMENSION OF A
C    N    I NUMBER OF ELEMENTS IN X = FIRST DIMENSION OF A
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING
C           (SUB)PROGRAM (L.GE.N)
C    A    I REAL ARRAY OF DIMENSION (L,K) WHERE K.GE.M
C    X    I GIVEN VECTOR, REAL VECTOR OF DIMENSION AT LEAST N
C    B    I RESULT, REAL VECTOR OF DIMENSION AT LEAST M
C
C    *******************************************************************
      SUBROUTINE MATTMV(M,N,L,A,X,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(L,M),X(N),B(M)
      DO 80 I=1,M
         T=0D0
         DO 40 J=1,N
  40        T=T+A(J,I)*X(J)
  80     B(I)=T
      END
