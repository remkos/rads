C    ****************************************************************** 
C                                                      35. LINMGS
C                                                                       
C    PURPOSE                                                            
C    TO SOLVE M LINEAR EQUATIONS IN N UNKNOWNS, WHERE M IS AT LEAST N
C    USING THE LEAST SQUARES METHOD AND ORTHOGONAL VECTORS
C                                                                       
C    USAGE                                                              
C    CALL LINMGS(A,X,B,Y,C,P,M,N,L,LA,LX,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I THE N GIVEN COLUMN VECTORS, DIMENSION (M,N)
C    X    O THE L RESULTANT VECTORS, DIMENSION (N,L)
C    B    I THE L GIVEN RIGHT-HAND SIDES, DIMENSION (M,L)
C    Y    O SQUARES OF THE LENGTHS OF PIVOTVECTORS
C    C      A MATRIX FOR WORKING SPACE OF DIMENSION AL LEAST (N,N)
C    P      A INTEGER VECTOR FOR WORKING SPACE OF DIMENSION AT LEAST N
C    M    I NUMBER OF ELEMENTS IN THE VECTORS OF A
C    N    I NUMBER OF VECTORS IN A
C    L    I NUMBER OF RIGHT HAND SIDES IN B
C    LA   I FIRST DIMENSION OF A AND B AS DECLARED IN THE CALLING
C           (SUB)PROGRAM. LA IS AT LEAST M
C    LX   I FIRST DIMENSION OF C AND X AS DECLARED IN THE CALLING
C           (SUB)PROGRAM. LX IS AT LEAST N
C    EP   I REAL TOLERANCE VALUE                                        
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT                         
C                              IF=1 PIVOT LESS EP
C                                                                       
C    METHOD                                                             
C    LEAST SQUARES AND ORTHOGONAL VECTORS
C                                                                       
C    *******************************************************************
      SUBROUTINE LINMGS(A,X,B,Y,C,P,M,N,L,LA,LX,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER P,Q
      DIMENSION A(LA,N),X(LX,N),B(LA,L),Y(N),C(LX,N),P(N)
      DATA Z,U/0.0D0,1.0D0/
      IF=0
      DO 20 K=1,N
      S=Z
      P(K)=K
      DO 10 I=1,M
  10  S=S+A(I,K)*A(I,K)
  20  Y(K)=S
      DO 120 K=1,N
      S=Z
      DO 24 J=K,N
      IF (S.LT.Y(J)) THEN
      S=Y(J)
      Q=J
      ENDIF
  24  CONTINUE
      IF (Q.GT.K) THEN
      P(K)=Q
      DO 30 I=1,M
      T=A(I,K)
      A(I,K)=A(I,Q)
  30  A(I,Q)=T
      DO 40 I=K-1,1,-1
      T=C(I,K)
      C(I,K)=C(I,Q)
  40  C(I,Q)=T
      T=Y(K)
      Y(K)=Y(Q)
      Y(Q)=T
      ENDIF
      IF (K.GT.1) THEN
      T=Z
      DO 50 I=1,M
  50  T=T+A(I,K)*A(I,K)
      ELSE
      T=S
      ENDIF
      Y(K)=T
      IF (T.LT.EP) THEN
      IF=1
      GOTO 200
      ENDIF
      DO 80 J=K+1,N
      S=Z
      DO 60 I=1,M
  60  S=S+A(I,J)*A(I,K)
      S=S/T
      C(K,J)=S
      DO 70 I=1,M
  70  A(I,J)=A(I,J)-A(I,K)*S
  80  Y(J)=Y(J)-S*S*T
      DO 110 J=1,L
      S=Z
      DO 90 I=1,M
  90  S=S+B(I,J)*A(I,K)
      S=S/T
      X(K,J)=S
      DO 100 I=1,M
 100  B(I,J)=B(I,J)-A(I,K)*S
 110  CONTINUE
 120  CONTINUE
      DO 140 K=1,L
      DO 130 I=N,1,-1
      S=Z
      DO 124 J=I+1,N
 124  S=S+C(I,J)*X(J,K)
 130  X(I,K)=X(I,K)-S
 140  CONTINUE
      DO 160 K=N-1,1,-1
      I=P(K)
      IF (I.EQ.K) GOTO 160
      DO 150 J=1,L
      T=X(I,J)
      X(I,J)=X(K,J)
 150  X(K,J)=T
 160  CONTINUE
 200  RETURN
      END
