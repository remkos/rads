C    ******************************************************************
C                                                      48. MATSSM
C
C    PURPOSE
C    TO MULTIPLY A MATRIX BY A MATRIX  AB=C,
C    WHERE A AND B ARE SYMMETRIC MATRICES IN VECTOR FORM
C
C    USAGE
C    CALL MATSSM(N,L,H,A,B,C)
C
C    DESCRIPTION OF PARAMETERS
C    N    I DIMENSION OF ALL MATRICES
C    L    I FIRST DIMENSION OF C AS DECLARED IN THE CALLING
C           (SUB)PROGRAM (L.GE.N)
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH SUBROUTINE MATSY1
C    A,B  I REAL VECTORS CONTAINING THE COLUMS OF THE UPPER PART OF
C           EACH SYMMETRIC MATRIX. THE DIMENSION IS AT LEAST N*(N+1)/2
C           A AND B CAN BE FORMED BY SUBROUTINE MATSYD
C    C    O RESULTANT MATRIX OF DIMENSION (L,N)
C
C    *******************************************************************
      SUBROUTINE MATSSM(N,L,H,A,B,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER H(N)
      DIMENSION A(*),B(*),C(L,N)
      DO 100 I=1,N
      DO 40 J=1,I
      T=0.0D0
      DO 10 K=1,J
  10  T=T+A(K+H(I))*B(K+H(J))
      DO 20 K=J+1,I
  20  T=T+A(K+H(I))*B(J+H(K))
      DO 30 K=I+1,N
  30  T=T+A(I+H(K))*B(J+H(K))
  40  C(I,J)=T
      DO 80 J=I+1,N
      T=0.0D0
      DO 50 K=1,I
  50  T=T+A(K+H(I))*B(K+H(J))
      DO 60 K=I+1,J
  60  T=T+A(I+H(K))*B(K+H(J))
      DO 70 K=J+1,N
  70  T=T+A(I+H(K))*B(J+H(K))
  80  C(I,J)=T
 100  CONTINUE
      END
