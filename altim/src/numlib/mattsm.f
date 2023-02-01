C    ******************************************************************
C                                                     46D. MATTSM
C
C    PURPOSE                      T
C    TO MULTIPLY THREE MATRICES  A BA=C
C    WHERE B AND C ARE SYMMETRIC MATRICES IN VECTOR FORM
C
C    USAGE
C    CALL MATTSM(N,M,L,H,V,A,B,C)
C
C    DESCRIPTION OF PARAMETERS
C    N    I DIMENSION OF B = FIRST DIMENSION OF A
C    M    I DIMENSION OF C = SECOND DIMENSION OF A
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING
C           (SUB)PROGRAM (L.GE.N)
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH SUBROUTINE MATSY1
C    V      WORKING SPACE OF DIMENSION AT LEAST N
C    A    I GIVEN MATRIX OF DIMENSION (L,M)
C    B    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX B. THE DIMENSION IS AT LEAST N*(N+1)/2
C           B CAN BE FORMED BY SUBROUTINE MATSYD
C    C    O RESULTANT REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX C. THE DIMENSION IS AT LEAST M*(M+1)/2
C           C CAN BE RETRANSFORMED BY SUBROUTINE MATSYR
C
C    *******************************************************************
      SUBROUTINE MATTSM(N,M,L,H,V,A,B,C)
      INTEGER H(N),S,R
      REAL*8 V(N),A(L,M),B(*),C(*),T
      K=0
      DO 80 I=1,M
         DO 40 R=1,N
            T=0
            DO 10 S=1,R
   10          T=T+A(S,I)*B(S+H(R))
            DO 20 S=R+1,N
   20          T=T+A(S,I)*B(R+H(S))
   40       V(R)=T
         DO 80 J=1,I
            K=K+1
            T=0
            DO 50 R=1,N
   50          T=T+V(R)*A(R,J)
   80       C(K)=T
      END
