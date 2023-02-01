C    ******************************************************************
C                                                      47. MATSMV
C
C    PURPOSE
C    TO MULTIPLY A VECTOR BY A MATRIX  AX=B,
C    WHERE A IS A SYMMETRIC MATRIX IN VECTOR FORM
C
C    USAGE
C    CALL MATSMV(N,H,A,X,B)
C
C    DESCRIPTION OF PARAMETERS
C    N    I NUMBER OF ELEMENTS IN X AND B
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH SUBROUTINE MATSY1
C    A    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX. THE DIMENSION IS AT LEAST N*(N+1)/2
C    X    I GIVEN VECTOR, REAL VECTOR OF DIMENSION AT LEAST N
C    B    O RESULT, REAL VECTOR OF DIMENSION AT LEAST M
C
C    *******************************************************************
      SUBROUTINE MATSMV(N,H,A,X,B)
      INTEGER H(N)
      REAL*8 A(*),X(N),B(N)
      DO 10 I=1,N
  10     B(I)=0D0
      DO 30 J=1,N
         DO 20 I=1,J
  20        B(I)=B(I)+A(I+H(J))*X(J)
         DO 30 I=J+1,N
  30        B(I)=B(I)+A(J+H(I))*X(J)
      END
