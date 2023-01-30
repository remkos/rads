C    ****************************************************************** 
C                                                      61. RUNEEN
C                                                                       
C    PURPOSE                                                            
C    TO INTEGRATE A FIRST ORDER ORDINARY DIFFERENTIAL EQUATION
C    Y'=F(X,Y), WHERE Y(X0)=Y0
C                                                                       
C    USAGE                                                              
C    CALL RUNEEN(X,Y,XN,N,F,C,D)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    X    I INITIAL VALUE X0 OF THE ARGUMENT X
C         O AT THE END X HAS THE VALUE XN
C    Y    I Y0, INITIAL VALUE Y0(X0) AT POINT X0
C         O FUNCTION VALUE AT XN
C    XN   I UPPER BOUND OF THE INTGRATION INTERVAL
C    N    I NUMBER OF INTGRATION STEPS, STEP LENGTH IS (XN-X)/N
C    F    I NAME OF DOUBLE PRECISION FUNCTION F(X,Y) TO BE INTEGRATED
C    C    O ARRAY C(N) CONTAINS THE ARGUMENT VALUES XI=X0+I*H
C    D    O ARRAY D(N) CONTAINS THE FUNCTION VALUES YI AT XI
C                                                                       
C    METHOD                                                             
C    FOURTH-ORDER RUNGE-KUTTA-GILL INTEGRATION PROCESS
C
C    *******************************************************************
      SUBROUTINE RUNEEN(X,Y,XN,N,F,C,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION K
      DIMENSION C(N),D(N)
      EXTERNAL F
      DATA P,A1,A2,A3,A4/0.333333333333333D0,0.585786437626905D0,3.41421
     13562373095D0,0.1213203435596426D0,-4.121320343559643D0/
      H=(XN-X)/N
      H2=H/2.0D0
      X0=X
      DO 20 I=1,N
      K=H2*F(X,Y)
      Y=Y+K
      Q=K
      K=H2*F(X+H2,Y)
      Y=Y+A1*(K-Q)
      Q=A1*K+A3*Q
      K=H2*F(X+H2,Y)
      Y=Y+A2*(K-Q)
      Q=A2*K+A4*Q
      K=H2*F(X+H,Y)
      Y=Y+P*(K-Q-Q)
      X=X0+I*H
      C(I)=X
      D(I)=Y
  20  CONTINUE
      RETURN
      END
