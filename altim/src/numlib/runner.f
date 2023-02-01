C    ****************************************************************** 
C                                                      62. RUNNER
C                                                                       
C    PURPOSE                                                            
C    TO INTEGRATE N-1 FIRST ORDER ORDINARY DIFFERENTIAL EQUATION
C     DXI/DX1 = FI(X1,X2,...XN+1), I=2(1)N
C     WITH INITIAL VALUES XI0 AT X10
C                                                                       
C    USAGE                                                              
C    CALL RUNNER(X,N,H,F,S,P,K)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    X    I INITIAL VALUES OF THE VARIABLES
C         O COMPUTED VALUES OF THESE VARIABLES AFTER ONE STEP H
C    N    I NUMBER OF EQUATIONS
C    H    I STEPLENGTH
C    F    I SUBROUTINE GIVING THE RIGHT HAND SIDE OF THE SYSTEM
C    S      REAL WORKING ARRAY OF DIMENSION (N)
C    P      REAL WORKING ARRAY OF DIMENSION (N)
C    K      REAL WORKING ARRAY OF DIMENSION (4,N)
C
C    REMARKS
C    SUBROUTINE F(G,X) YIELDS IN ARRAY G THE N FUNCTION VALUES
C    FI(X) AT X10, I=2(1)N
C                                                                       
C    METHOD                                                             
C    FOURTH-ORDER RUNGE-KUTTA INTEGRATION PROCESS
C
C    *******************************************************************
      SUBROUTINE RUNNER(X,N,H,F,S,P,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION K
      DIMENSION X(N),S(N),P(N),K(4,N)
      EXTERNAL F
      DO 10 I=1,N
  10  P(I)=X(I)
      I=1
      S(1)=1.D0
  20  CALL F(S,P)
      DO 40 J=1,N
      K(I,J)=H*S(J)
      IF (I.LT.3) THEN
      P(J)=X(J)+K(I,J)/2.D0
      ELSE
      P(J)=X(J)+K(I,J)
      ENDIF
  40  CONTINUE
      IF (I.LT.4) THEN
      I=I+1
      GOTO 20
      ENDIF
      DO 60 J=1,N
      P(J)=X(J)+(K(1,J)+2.D0*(K(2,J)+K(3,J))+K(4,J))/6.D0
      X(J)=P(J)
  60  CONTINUE
      RETURN
      END
