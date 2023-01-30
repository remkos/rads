C    ****************************************************************** 
C                                                      18. INTEXP
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE INTEGRAL OF A FUNCTION F*DEXP(-X)
C    ON INTERVAL (X0,INF)
C                                                                       
C    USAGE                                                              
C    CALL INTEXP(F,X0,HI,EP,T,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    F    I NAME OF DOUBLE PRECISION FUNCTION TO BE INTEGRATED
C    X0   I LOWER BOUND OF THE INTEGRAL
C    HI   I INITIAL INCREMENT OF THE ARGUMENT VALUES
C    EP   I ABSOLUTE TOLERANCE VALUE FOR THE SUBINTERVALS
C    T    O THE RESULTANT INTEGRAL VALUE
C    IF   O INTEGER ERRORCODE IF = 0 NORMAL EXIT
C                             IF = 1 N EXEEDS NMAX
C                                                                       
C    METHOD                                                             
C    GAUSS-LAGUERRE QUADRATURE FORMULAS
C
C    REMARKS
C    INTGAU USES FUNCTION INTGAU, INTEX1 AND INTEX2
C                                                                       
C    *******************************************************************
      SUBROUTINE INTEXP(F,X0,HI,EP,T,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION INTEX1,INTEX2
      DIMENSION A(10),H(10)
      EXTERNAL F,INTEX2
      DATA A/0.13779347054D0,0.7294545495D0,1.80834290274D0,3.4014336978
     16D0,5.55249614006D0,8.33015274676D0,1.18437858379D1,1.62792578314D
     21,2.1996585812D1,2.99206970123D1/
      DATA H/0.308441115765D0,0.401119929155D0,0.218068287612D0,0.620874
     1560987D-1,0.950151697518D-2,0.753008388588D-3,0.28259233496D-4,0.4
     224931398496D-6,0.183956482398D-8,0.991182721961D-12/
      B=X0
      Q=HI
      N=1000
      T=0.0D0
      IF=0
      IF (Q.LT.0001) Q=1.0D0
      P0=INTEX1(F,B,A,H)
      DO 60 I=1,100
      S=B+I*Q
      GOTO 20
  10  EP=2.D0*EP
  20  IP=0
      CALL INTGAU(S-Q,S,Q,EP,F,INTEX2,K,N,D,IP)
      IF (IP.EQ.1) GOTO 10
      T=T+D
      P1=INTEX1(F,S,A,H)
      IF (DABS(P1+D-P0).LT.EP*DABS(T+P1)) GOTO 80
      P0=P1
  60  CONTINUE
      IF=1
  80  T=T+P1
 100  FORMAT(A,I6,3F24.16)
      RETURN
      END
      DOUBLE PRECISION FUNCTION INTEX1(F,B,A,H)
      DOUBLE PRECISION F,B,A,H,T
      DIMENSION A(10),H(10)
      EXTERNAL F
      T=0.D0
      DO 10 I=1,10
  10  T=T+H(I)*F(A(I)+B)
      IF (B.GT.1.D-12) T=T*DEXP(-B)
      INTEX1=T
      RETURN
      END
      DOUBLE PRECISION FUNCTION INTEX2(X,H,F)
      DOUBLE PRECISION X,Y,H,R,T,F,YE,RE
      EXTERNAL F
      Y=X+H/2.D0
      YE=DEXP(-Y)
      R=0.453089922969332D0*H
      RE=DEXP(-R)
      T=0.1184634425280945D0*(F(Y+R)*RE+F(Y-R)/RE)
      R=0.2692346550528415D0*H
      RE=DEXP(-R)
      T=T+0.239314335249683D0*(F(Y+R)*RE+F(Y-R)/RE)
      INTEX2=H*(T+0.2844444444444444D0*F(Y))*YE
      RETURN
      END
