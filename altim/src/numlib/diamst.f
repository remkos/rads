C    ******************************************************************
C                                                     46C. DIAMST
C
C    PURPOSE                                            T
C    TO COMPUTE THE DIAGONAL OF A SYMMETRIC MATRIX C=ABA
C    WHERE B IS A SYMMETRIC MATRIX IN VECTOR FORM
C
C    USAGE
C    CALL DIAMST(N,M,L,H,A,B,C)
C
C    DESCRIPTION OF PARAMETERS
C    N    I DIMENSION OF B = SECOND DIMENSION OF A
C    M    I DIMENSION OF C = FIRST DIMENSION OF A
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING
C           (SUB)PROGRAM (L.GE.M)
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH SUBROUTINE MATSY1
C    A    I GIVEN MATRIX OF DIMENSION (L,N)
C    B    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX B. THE DIMENSION IS AT LEAST N*(N+1)/2
C           B CAN BE FORMED BY SUBROUTINE MATSYD
C    C    O RESULTANT REAL VECTOR CONTAINING THE DIAGONAL OF
C           THE SYMMETRIC MATRIX C.
C
C    *******************************************************************
      SUBROUTINE DIAMST(N,M,L,H,A,B,C)
      INTEGER H(N),S,R
      REAL*8 A(L,N),B(*),C(M),T,TT
      K=0
      DO 80 I=1,M
         TT=0
         DO 40 R=1,N
            T=0
            DO 10 S=1,R
   10          T=T+A(I,S)*B(S+H(R))
            DO 20 S=R+1,N
   20          T=T+A(I,S)*B(R+H(S))
   40       TT=TT+T*A(I,R)
   80    C(I)=TT
      END
