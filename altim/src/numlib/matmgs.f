C    ****************************************************************** 
C                                                      43. MATMGS
C                                                                       
C    PURPOSE                                                            
C    TO ORTHOGONALIZE A SYSTEM OF VECTORS BY
C    MEANS OF THE MODIFIED GRAM-SCHMIDT METHOD
C                                                                       
C    USAGE                                                              
C    CALL MATMGS(A,C,D,M,N,LA,LC,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I N COLUMN VECTORS DIMENSION (M,N)
C         O THE N RESULTANT ORTHONORMAL VECTORS
C    C    O A TRIANGULAR MATRIX OF COEFFICIENTS OF DIMENSION
C           (N,N) OF WHICH THE UPPER TRIANGEL IS USED
C    D    O PRODUCT OF THE MAIN DIAGONAL ELEMENTS OF C
C           IF M=N THEN D IS THE DETERMINANT OF A
C    M    I NUMBER OF ELEMENTS IN THE VECTORS OF A
C    N    I NUMBER OF VECTORS IN A
C    LA   I FIRST DIMENSION OF A AS DECLARED IN THE CALLING
C           (SUB)PROGRAM. LA IS AT LEAST M
C    LC   I FIRST DIMENSION OF C AS DECLARED IN THE CALLING
C           (SUB)PROGRAM. LC IS AT LEAST N
C    EP   I REAL TOLERANCE VALUE                                        
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT                         
C                              IF=1 PIVOT LESS EP
C                                                                       
C    METHOD                                                             
C    MODIFIED GRAM-SCHMIDT METHOD
C                                                                       
C    *******************************************************************
      SUBROUTINE MATMGS(A,C,D,M,N,LA,LC,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LA,M),C(LC,N)
      DATA Z,U/0.0D0,1.0D0/
      IF=0
      D=U
      DO 50 K=1,N
      T=Z
      DO 10 I=1,M
  10  T=T+A(I,K)*A(I,K)
      C(K,K)=DSQRT(T)
      IF (C(K,K).LT.EP) THEN
      D=C(K,K)
      IF=1
      GOTO 100
      ENDIF
      D=D*C(K,K)
      DO 40 J=K+1,N
      S=Z
      DO 20 I=1,M
  20  S=S+A(I,J)*A(I,K)
      S=S/T
      C(K,J)=S
      DO 30 I=1,M
  30  A(I,J)=A(I,J)-A(I,K)*S
  40  CONTINUE
  50  CONTINUE
      DO 80 J=1,N
      T=C(J,J)
      S=U/T
      DO 60 I=1,M
  60  A(I,J)=A(I,J)*S
      DO 70 I=J+1,N
  70  C(J,I)=C(J,I)*T
  80  CONTINUE
 100  RETURN
      END
