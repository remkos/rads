C    ****************************************************************** 
C                                                      12. ERRORF
C                                                                       
C    PURPOSE                                                            
C    TO INTEGRATE THE ERROR FUNCTION ERF(X):
C     INTEGRAL (DT:0,X) (2/SQRT(PI) * E**(-T**2))
C    AND THE COMPLEMENTARY ERROR FUNCTION ERFC(X) = 1-ERF(X)
C                                                                       
C    USAGE                                                              
C    CALL ERRORF(X,ERF,ERFC,L)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    X    I ARGUMENT VALUE
C    ERF  O THE VALUE OF ERF(X)
C    ERFC O THE VALUE OF ERFC(X)
C           THE VALUE ERFC(X)-2 IF X<0
C    L    O L=1  IF X>=0
C           L=-1 IF X<0
C
C    REMARKS
C    IF X>5 THEN 1-ERF(X)<10**-16. THEN ERFC(X) CAN BE USED.
C    IF X>=0 THEN ERFC(X) = 1-ERF(X)
C    IF X<0  THEN ERFC(X) = (-1)-ERF(X) = -(1+ERF(X))
C    ERRORF USES SUBROUTINE ERROR1(N,X,A,Y)
C
C    METHOD                                                             
C    REFER TO : C.W. CLENSHAW, MATHEMATICAL TABLES, VOL 5 PAGE 28
C
C    *******************************************************************
      SUBROUTINE ERRORF(X,ERF,ERFC,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(55),B(27)
      DATA A/1.9436518276114522D0,0.D0,-1.3816314200197992D0,0.D0,0.6473
     1164048545842D0,0.D0,-.30593102442203560D0,0.D0,0.1386797472020301D
     20,0.D0,-0.59247456591259D-1,0.D0,0.236917518249282D-1,0.D0,-0.8847
     33626352405D-2,0.D0,.30856617113609D-2,0.D0,-.1006386351238D-2,0.D0
     4,.3075463288431D-3,0.D0,-0.882619837554D-4,0.D0,.238450961661D-4,0
     5.D0,-0.60791002851D-5,0.D0,0.14659721734D-5,0.D0,-0.3351599343D-6,
     60.D0,0.728057954D-7,0.D0,-0.150579118D-7,0.D0,0.29709474D-8,0.D0,-
     70.5602127D-9,0.D0,0.1011316D-9,0.D0,-0.175065D-10,0.D0,0.29104D-11
     8,0.D0,-0.4653D-12,0.D0,0.716D-13,0.D0,-0.106D-13,0.D0,0.15D-14,0.D
     90,-0.2D-15/
      DATA B/0.9853526361287725,0.D0,-0.14339740271775D-1,0.D0,0.2973616
     1922026D-3,0.D0,-0.98035160434D-5,0.D0,0.4331334203D-6,0.D0,-0.2362
     215003D-7,0.D0,0.15154968D-8,0.D0,-0.1108494D-9,0.D0,0.90426D-11,0.
     3D0,-0.8095D-12,0.D0,0.785D-13,0.D0,-0.82D-14,0.D0,0.1D-14,0.D0,-0.
     41D-15/
      L=-1
      IF (X.LE.0.D0) THEN
      X=-X
      ELSE
      L=1
      ENDIF
      IF (X.LE.4.D0) THEN
      S=X/4.D0
      CALL ERROR1(55,S,A,Y)
      ERF=S*Y*L
      ERFC=1-ERF
      ELSE
      S=4.0D0/X
      CALL ERROR1(27,S,B,Y)
      ERFC=0.5641895835477563D0*L*DEXP(-X*X)*Y/X
      ERF=1-ERFC
      ENDIF
      RETURN
      END
      SUBROUTINE ERROR1(N,X,A,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N)
      U=X+X
      PB=0.D0
      PA=0.D0
      DO 20 I=N,2,-1
      PC=U*PB+A(I)-PA
      PA=PB
  20  PB=PC
      Y=PB*X-PA+A(1)
      RETURN
      END
