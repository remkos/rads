C    ****************************************************************** 
C                                                      19. INTGAU
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE A DEFINITE INTEGRAL OF A FUNCTION F
C    ON INTERVAL (X0,XN)
C                                                                       
C    USAGE                                                              
C    CALL INTGAU(X0,XN,HI,EP,F,INT,N,NMAX,S,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    X0   I LOWER BOUND OF THE INTEGRAL
C    XN   I UPPER BOUND OF THE INTEGRAL
C    HI   I INITIAL INCREMENT OF THE ARGUMENT VALUES
C    EP   I ABSOLUTE TOLERANCE VALUE FOR THE SUBINTERVALS
C    F    I NAME OF DOUBLE PRECISION FUNCTION TO BE INTEGRATED
C    INT  I FUNCTION TO EVALUATE QUADRATURE FORMULA
C           INTGA1 MUST BE PUT HERE
C    N    O NUMBER OF SUBINTERVALS USED
C    NMAX I MAXIMUM NUMBER OF SUBINTERVALS
C    S    O THE RESULTANT INTEGRAL VALUE
C    IF   O INTEGER ERRORCODE IF = 0 NORMAL EXIT
C                             IF = 1 N EXEEDS NMAX
C                                                                       
C    METHOD                                                             
C    5-POINT GAUSSIAN-LEGENDRE QUADRATURE FORMULA WITH SUBINTERVALS
C
C    REMARKS
C    INTGAU USES FUNCTION INTGA1
C                                                                       
C    *******************************************************************
      SUBROUTINE INTGAU(X0,XN,HI,EP,F,INT,N,NMAX,S,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION INT
      LOGICAL B,B1
      DIMENSION FF(7),HH(7)
      EXTERNAL F,INT
      DATA Z/0.0D0/
      H=HI
      XA=X0
      XB=XN
      NM=NMAX
      T=XB-XA
      B1=(T.LT.Z)
      IF (B1) THEN
      R=XA
      XA=XB
      XB=R
      T=-T
      ENDIF
      IF (H.LT.1.D-10) H=T
      S=Z
      B=.FALSE.
      HH(1)=H
      H1=H
      N=0
      K=1
      IF=0
      DO 200 I=1,NM
      IF (K.EQ.1) THEN
      IF (B) GOTO 210
      K=2
      IF (XA+1.001D0*H1.GT.XB) THEN
      H1=XB-XA
      B=.TRUE.
      ENDIF
      HH(2)=H1
      HH(1)=H1+H1
      FF(2)=INT(XA,H1,F)
      N=N+1
      ENDIF
      H2=H1/2.D0
      T1=INT(XA,H2,F)
      T2=INT(XA+H2,H2,F)
      N=N+2
      T=T1+T2
      IF (DABS(T-FF(K)).LT.EP) THEN
      S=S+T
      XA=XA+H1
      K=K-1
      H1=HH(K)
      ELSE
      H1=H2
      FF(K)=T2
      HH(K)=H1
      K=K+1
      IF (K.GT.7) THEN
      B=.FALSE.
      K=K-1
      DO 100 J=2,K
      FF(J-1)=FF(J)
  100 HH(J-1)=HH(J)
      ENDIF
      FF(K)=T1
      HH(K)=H1
      ENDIF
  200 CONTINUE
      IF=1
  210 IF (B1) S=-S
      RETURN
      END
