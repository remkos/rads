C    ****************************************************************** 
C                                                      22. INVERF
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE INVERSE OF THE ERROR FUNCTION ERF(X):
C     INTEGRAL (DT:0,X) (2/SQRT(PI) * E**(-T**2))
C    OR THE INVERSE OF
C     THE COMPLEMENTARY ERROR FUNCTION ERFC(X) = 1-ERF(X)
C                                                                       
C    USAGE                                                              
C    CALL INVERF(IS,X,Y,IC)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    IS   I IS=0 INVERF(X)  IS COMPUTED
C           IS=1 INVERFC(X) IS COMPUTED
C    X    I ARGUMENT VALUE
C    Y    O FUNCTION VALUE IN ACCORDANCE WITH OPTIONS IS
C    IC   O RETURN CODE : DETERMINED BY THE VALUES OF X AND IS
C           SEE FOLLOWING TABLE, WHERE X>0 IS ASSUMED
C
C    ------------------------------------------------------------
C    IC |             IS=0            |           IS=1
C    ------------------------------------------------------------
C    1    X<=0.8                         0.2<=X<=1
C    2    0.8<X<=1-25*10**-4             25*10**-4<=X<0.2
C    3    1-25*10**-4<X<=1-5*10**-16     5*10**-16<=X<25*10**-4
C    4                                   10**-74<=X<5*10**-16
C    5    1-5*10**-16<X<=1, Y=6
C    6                                   X<10**-74, Y=13
C    7    X>1, Y=0                       X>1, Y=0
C    ------------------------------------------------------------
C    IC>0 IF X>=0 AND IC<0 IF X<0
C
C    REMARKS
C
C    IF X>1.5*10**-16 THEN INERF(X)=1. IN THIS CASE THE INVERSE
C    COMPLEMENTARY ERROR FUNCTION CAN STILL BE USED
C    FOR IS=1 AND X>0 Y=INVERF(1-X)=INVERFC(X)
C    FOR IS=1 AND X<0 Y=INVERF(-1-X)
C    INVERF USES SUBROUTINE ERROR1
C
C    METHOD
C    A.J. STRECOK, ON THE CALCULATION OF THE INVERSE OF THE
C    ERROR FUNCTION.
C    *******************************************************************
      SUBROUTINE INVERF(IS,X,Y,IC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(24),B(18),C(25),D(14)
      DATA A/0.9928853766189408D0,0.1204675161431045D0,0.160781993420999
     1D-1,0.26867044371623D-2,0.4996347302357D-3,0.988982185991D-4,0.203
     291812764D-4,0.43272716177D-5,0.9380814129D-6,0.2067347209D-6,0.461
     3596991D-7,0.104166797D-7,0.2371501D-8,0.5439284D-9,.125549D-9,0.29
     41382D-10,0.67949D-11,0.15912D-11,0.374D-12,0.882D-13,0.209D-13,0.4
     59D-14,0.12D-14,0.3D-15/
      DATA B/0.9121588034175538D0,-0.162662818676637D-1,0.4335564729494D
     1-3,0.2144385700745D-3,0.26257510758D-5,-0.30210910501D-5,-0.124060
     2619D-7,0.624066093D-7,-0.5401248D-9,-0.14232079D-8,0.34384D-10,0.3
     335849D-10,-0.14584D-11,-0.8102D-12,0.525D-13,0.197D-13,-0.17D-14,-
     40.5D-15/
      DATA C/0.9566797090204925D0,-0.231070043090649D-1,-0.4374236097508
     14D-2,-0.5765043226512D-3,-0.109610223071D-4,0.251085470246D-4,0.10
     25623360679D-4,0.275441233D-5,0.4324844983D-6,-0.205303367D-7,-0.43
     38915367D-7,-0.176840095D-7,-0.3991289D-8,-0.1869324D-9,0.2729227D-
     49,0.1328172D-9,0.318342D-10,0.16701D-11,-0.20365D-11,-0.9648D-12,-
     50.2196D-12,-0.96D-14,0.137D-13,0.63D-14,0.15D-14/
      DATA D/0.9885750640661893D0,0.108577051845995D-1,-0.17511651027628
     1D-2,0.211969932066D-4,0.156648714042D-4,-0.5190416869D-6,-0.371357
     2879D-7,0.12174309D-8,-0.1768116D-9,-0.119372D-10,0.3803D-12,-0.66D
     3-13,-0.88D-14,-0.4D-15/
      XX=X
      L=-1
      IF (X.LT.0.D0) THEN
      X=-X
      ELSE
      L=1
      ENDIF
      IF (X.GT.1.D0) THEN
      Y=0.D0
      IC=7
      GOTO 200
      ENDIF
      IF (IS.EQ.1) THEN
      U=X
      ELSE
      U=1.D0-X
      ENDIF
      IF (U.GE.0.2D0) GOTO 10
      IF (U.GE.1.D-74) BETA=DSQRT(-DLOG(U*(2.D0-U)))
      IF (U.GE.25.D-4) GOTO 20
      IF (U.GE.5.D-16) GOTO 30
      IF (IS.EQ.0) GOTO 50
      IF (U.GE.1.D-74) GOTO 40
      Y=13.D0
      IC=6
      GOTO 200
  10  U=1.D0-U
      W=U*U/0.32D0-1.D0
      CALL ERROR1(24,W,A,Y)
      Y=Y*U
      IC=1
      GOTO 200
  20  IC=2
      W=-1.548813042373262D0*BETA+2.565490123147816D0
      CALL ERROR1(18,W,B,Y)
      GOTO 100
  30  IC=3
      W=-0.559457631329832D0*BETA+2.287915716263358D0
      CALL ERROR1(25,W,C,Y)
      GOTO 100
  40  IC=4
      W=-9.199992358830151D0/DSQRT(BETA)+2.794990820124599D0
      CALL ERROR1(14,W,D,Y)
      GOTO 100
  50  Y=6.D0
      IC=5
      GOTO 200
 100  Y=Y*BETA
 200  Y=Y*L
      IC=IC*L
      X=XX
      RETURN
      END
