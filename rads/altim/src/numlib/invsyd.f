C    ****************************************************************** 
C                                                      24. INVSYD
C                                                                       
C    PURPOSE                                                            
C    TO INVERT A MATRIX A
C    WHERE A IS A SYMMETRIC MATRIX IN VECTOR FORM
C                                                                       
C    USAGE                                                              
C    CALL INVSYD(N,A,H,P,Q,R,D,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I DIMENSION OF A
C    A    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX A. THE DIMENSION IS AT LEAST N*(N+1)/2
C           A CAN BE FORMED BY SUBROUTINE MATSYD
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH MATSY1
C    P,Q    REAL VECTORS, WORKING SPACE OF DIMENSION AT LEAST N
C    R      LOGICAL VECTOR WORKING SPACE OF DIMENSION AT LEAST N
C    D    O DETERMINANT OF A                                            
C    EP   I REAL TOLERANCE VALUE                                        
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT                         
C                              IF=1 PIVOT LESS EP
C                                                                       
C    METHOD                                                             
C    GAUSS ELIMINATION WITH PIVOTTING ON THE MAIN DIAGONAL
C                                                                       
C    *******************************************************************
      SUBROUTINE INVSYD(N,A,H,P,Q,R,D,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER H
      LOGICAL R
      DIMENSION A(*),H(N),P(N),Q(N),R(N)
      DATA Z,U/0.0D0,1.0D0/
      IF=0
      D=U
      DO 2 I=1,N
   2  R(I)=.TRUE.
      DO 30 I=1,N
      Y=Z
      DO 4 J=1,N
      T=DABS(A(J+H(J)))
      IF ((Y.LT.T) .AND. (R(J))) THEN
      K=J
      Y=T
      ENDIF
   4  CONTINUE
      IF (Y.GE.EP) GOTO 10
      IF=1
      GOTO 100
  10  R(K)=.FALSE.
      M=H(K)
      Q(K)=U/A(K+M)
      D=D*A(K+M)
      P(K)=U
      A(K+M)=Z
      DO 15 J=K-1,1,-1
      T=A(J+M)
      P(J)=T
      Q(J)=Q(K)*T
      IF (R(J)) Q(J)=-Q(J)
  15  A(J+M)=Z
      DO 20 J=K+1,N
      M=H(J)
      T=A(K+M)
      P(J)=-T
      IF (R(J)) P(J)=T
      Q(J)=-Q(K)*T
  20  A(K+M)=Z
      DO 24 K=1,N
      M=H(K)
      DO 22 J=1,K
  22  A(J+M)=A(J+M)+P(J)*Q(K)
  24  CONTINUE
  30  CONTINUE
 100  RETURN
      END
