C    ****************************************************************** 
C                                                      23. INVGAU
C                                                                       
C    PURPOSE                                                            
C    TO INVERT A SQUARE MATRIX A
C                                                                       
C    USAGE                                                              
C    CALL INVGAU(N,L,A,P,Q,B,C,D,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I NUMBER OF EQUATIONS                                         
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING             
C           (SUB)PROGRAM (L.GE.N)                                       
C    A    I REAL ARRAY OF DIMENSION (L,K) WHERE K.GE.N
C         O THE RESULTANT INVERSE MATRIX
C    P    I INTEGER WORKING SPACE OF DIMENSION AT LEAST N
C    Q    I AS P
C    B    I REAL WORKING SPACE OF DIMENSION AT LEAST N
C    C    I AS B
C    D    O DETERMINANT OF A                                            
C    EP   I REAL TOLERANCE VALUE                                        
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT                         
C                              IF=1 PIVOT LESS EP
C                                                                       
C    METHOD                                                             
C    GAUSS ELIMINATION WITH COLUMN PIVOTTING                            
C                                                                       
C    *******************************************************************
      SUBROUTINE INVGAU(N,L,A,P,Q,B,C,D,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER P,Q
      DIMENSION A(L,N),P(N),Q(N),B(N),C(N)
      DATA Z,U/0.0D0,1.0D0/
      IF=0
      D=U
      DO 19 K=1,N
      PIV=Z
      DO 4 I=K,N
      DO 4 J=K,N
      T=DABS(A(I,J))
      IF (PIV.GE.T) GOTO 4
      P(K)=I
      Q(K)=J
      PIV=T
   4  CONTINUE
      IF (PIV.GE.EP) GOTO 6
      IF=1
      GOTO 100
   6  I=P(K)
      IF (K.EQ.I) GOTO 10
      DO 8 J=1,N
      T=A(I,J)
      A(I,J)=A(K,J)
   8  A(K,J)=T
      D=-D
  10  J=Q(K)
      IF (K.EQ.J) GOTO 14
      DO 12 I=1,N
      T=A(I,J)
      A(I,J)=A(I,K)
  12  A(I,K)=T
      D=-D
  14  D=D*A(K,K)
      R=U/A(K,K)
      DO 16 J=1,N
      IF (J.EQ.K) THEN
      B(J)=R
      C(J)=U
      ELSE
      B(J)=-A(K,J)*R
      C(J)=A(J,K)
      END IF
      A(K,J)=Z
  16  A(J,K)=Z
      DO 18 J=1,N
      DO 18 I=1,N
      A(I,J)=A(I,J)+C(I)*B(J)
  18  CONTINUE
  19  CONTINUE
      DO 36 K=N,1,-1
      J=P(K)
      IF (K.EQ.J) GOTO 22
      DO 20 I=1,N
      T=A(I,J)
      A(I,J)=A(I,K)
  20  A(I,K)=T
  22	I=Q(K)
      IF (K.EQ.I) GOTO 36
      DO 24 J=1,N
      T=A(I,J)
      A(I,J)=A(K,J)
  24  A(K,J)=T
  36  CONTINUE
 100  RETURN
      END
