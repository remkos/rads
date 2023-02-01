C    ****************************************************************** 
C                                                      63. RUNNOR
C                                                                       
C    PURPOSE                                                            
C    TO INTEGRATE 1 STEP OF A NTH ORDER ORDINARY DIFFERENTIAL EQUATION
C     Y(N) = F(X,Y,Y(1),......Y(N-1)).
C     WITH INITIAL VALUES FOR X,Y,Y(1),...,Y(N-1)
C                                                                       
C    USAGE                                                              
C    CALL RUNNOR(F,X,N1,H,S,P,K)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    F    I FUNCTION F(X) GIVING THE DIFFERENTIAL EQUATION
C    X    I ARRAY X(N1) CONTAINS THE INITIAL VALUES OF
C           X,Y,Y(1),......,Y(N-1)
C         O COMPUTED VALUES OF THESE VARIABLES AFTER ONE STEP H
C    N1   I THE DIMENSION OF ARRAY X, N1=N+1
C    H    I STEPLENGTH (IN X)
C    S      REAL WORKING ARRAY OF DIMENSION (N1)
C    P      REAL WORKING ARRAY OF DIMENSION (N1)
C    K      REAL WORKING ARRAY OF DIMENSION (4,N1)
C
C    REMARKS
C    FUNCTION F(X) EVALUATES THE EXPRESSION FOR OBTAINING Y(N),
C    THE NTH DERIVATIVE OF Y
C    TO TEST ACCURACY OF THE RESULT THE PROCEDURE MAY BE REPEATED
C    WITH A STEPLENGTH = 0.5*H.
C                                                                       
C    METHOD                                                             
C    FOURTH-ORDER RUNGE-KUTTA INTEGRATION PROCESS
C
C    *******************************************************************
      SUBROUTINE RUNNOR(F,X,N1,H,S,P,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION K
      DIMENSION X(N1),S(N1),P(N1),K(4,N1)
      EXTERNAL F
      N=N1-1
      DO 10 I=1,N1
  10  P(I)=X(I)
      S(1)=1.D0
      DO 20 I=2,N
  20  S(I)=X(I+1)
      DO 60 I=1,4
      DO 30 J=1,N
  30  K(I,J)=H*S(J)
      K(I,N1)=H*F(X)
      DO 40 J=1,N1
  40  X(J)=P(J)+0.5D0*K(I,J)*((I/3)+1.D0)
      DO 50 J=2,N
  50  S(J)=X(J+1)
  60  CONTINUE
      DO 70 J=1,N1
  70  X(J)=P(J)+(K(1,J)+2.D0*(K(2,J)+K(3,J))+K(4,J))/6.D0
      RETURN
      END
