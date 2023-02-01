C    ****************************************************************** 
C                                                       1. BESSEL
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE FOLLOWING BESSEL FUNCTIONS:
C       JN(X), YN(X), IN(X) OR KN(X) FOR N=0 OR 1
C                                                                       
C    USAGE                                                              
C    VAR = BESSEL(K,N,X,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    K    I K FIXES THE NATURE OF THE BESSEL FUNCTIONS
C           K=1   JN(X)
C           K=2   YN(X)
C           K=3   IN(X)
C           K=4   KN(X)
C    N    I THE ORDER OF THE BESSEL FUNCTION (0 OR 1)
C    X    I ARGUMENT FOR THE BESSEL FUNCTION
C    IF   O ERROR CODE IF=0 NORMAL EXIT
C                      IF=1 (X < 0)     AND  (K = 1)
C                           (X <= 0)    AND  (K = 2)
C                           (X < -3.75) AND  (K = 3)
C                           (X < =0)    AND  (K = 4)
C
C    REMARKS
C    THE MAXIMUM ERROR IS OF ORDER 1.D-7
C    IF IF=1, BESSEL TAKES THE VALUE 1.D70
C    BESSEL USES THE SUBROUTINES BESSE1, BESSE2 AND BESSE3
C
C    METHOD                                                             
C    REFER TO : C.W. CLENSHAW, MATHEMATICAL TABLES, VOL 5 TABLE 8
C
C    *******************************************************************
      FUNCTION BESSEL(K,N,X,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA Z,U/0.0D0,1.0D0/
      IF=1
      HA=1.D70
      IF ((K.EQ.1) .AND. (X.LT.Z)) GOTO 200
      IF ((K.EQ.2) .AND. (X.LE.Z)) GOTO 200
      IF ((K.EQ.3) .AND. (X.LT.-3.75D0)) GOTO 200
      IF ((K.EQ.4) .AND. (X.LE.Z)) GOTO 200
      IF=0
      IF (K.LE.2) THEN
      HA=BESSE1(K,N,X)
      GOTO 200
      ELSE
      HA=BESSE2(K,N,X)
      ENDIF
 200  BESSEL=HA
      RETURN
      END
      FUNCTION BESSE1(K,N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(9),B(6)
      DATA Z,U/0.0D0,1.0D0/
      IF (N.EQ.0) THEN
      IF (X.LE.4.D0) THEN
      T=0.625D-1*X*X
      A(1)=-0.5014415D-3
      A(2)=0.76771853D-2
      A(3)=-0.709253492D-1
      A(4)=0.4443584263D0
      A(5)=-1.7777560599D0
      A(6)=3.9999973021D0
      A(7)=-3.9999998721D0
      A(8)=U
      CALL BESSE3(8,A,T,HELP)
      IF (K.EQ.1) GOTO 200
      A(1)=0.891322D-4
      A(2)=-0.13508487D-2
      A(3)=0.148999271D-1
      A(4)=-0.1214187561D0
      A(5)=0.6694321484D0
      A(6)=-2.2331102234D0
      A(7)=3.6911388793D0
      A(8)=-1.6911374142D0
      A(9)=-0.5772156649D0
      CALL BESSE3(9,A,T,R)
      HELP=0.636619772*(HELP*DLOG(0.5D0*X)-R)
      GOTO 200
      ENDIF
      T=16.D0/(X*X)
      A(1)=-0.37043D-5
      A(2)=0.173565D-4
      A(3)=-0.487613D-4
      A(4)=0.17343D-3
      A(5)=-0.1753062D-2
      A(6)=0.3989422793D0
      B(1)=0.32312D-5
      B(2)=-0.142078D-4
      B(3)=0.342468D-4
      B(4)=-0.869791D-4
      B(5)=0.4564324D-3
      B(6)=-0.124669441D-1
      CALL BESSE3(6,A,T,P)
      CALL BESSE3(6,B,T,R)
      Q=4.D0/X*R
      Y=X-0.78539816434D0
      IF (K.EQ.1) THEN
      HELP=2.D0*DSQRT(U/X)*(P*DCOS(Y)-Q*DSIN(Y))
      GOTO 200
      ENDIF
      HELP=2.D0*DSQRT(U/X)*(P*DSIN(Y)+Q*DCOS(Y))
      GOTO 200
      ENDIF
      IF (N.EQ.1) THEN
      IF (X.LE.4.D0)THEN
      T=X*X*0.625D-1
      A(1)=-0.1289769D-3
      A(2)=0.22069155D-2
      A(3)=-0.236616773D-1
      A(4)=0.1777582922D0
      A(5)=-0.8888839649D0
      A(6)=2.6666660544D0
      A(7)=-3.999999971D0
      A(8)=1.9999999998D0
      CALL BESSE3(8,A,T,R)
      HELP=X/4.D0*R
      IF (K.EQ.1) GOTO 200
      A(1)=-0.10266368D-2
      A(2)=0.169921876D-1
      A(3)=-0.169108172D0
      A(4)=1.1418033012D0
      A(5)=-4.9105291148D0
      A(6)=1.16207891416D1
      A(7)=-1.07645472724D1
      A(8)=-0.6177253972D0
      A(9)=1.0000000004D0
      CALL BESSE3(9,A,T,R)
      HELP=0.636619772D0*(HELP*DLOG(0.5D0*X)-R/X)
      GOTO 200
      ENDIF
      T=16.D0/(X*X)
      A(1)=0.42414D-5
      A(2)=-0.20092D-4
      A(3)=0.580759D-4
      A(4)=-0.223203D-3
      A(5)=0.29218256D-2
      A(6)=0.3989422819D0
      B(1)=-0.36594D-5
      B(2)=0.1622D-4
      B(3)=-0.398708D-4
      B(4)=0.1064741D-3
      B(5)=-0.63904D-3
      B(6)=0.374008364D-1
      CALL BESSE3(6,A,T,P)
      CALL BESSE3(6,B,T,R)
      Q=4.D0/X*R
      Y=X-2.3561944902D0
      IF (K.EQ.1) THEN
      HELP=2.D0*DSQRT(U/X)*(P*DCOS(Y)-Q*DSIN(Y))
      GOTO 200
      ENDIF
      HELP=2.D0*DSQRT(U/X)*(P*DSIN(Y)+Q*DCOS(Y))
      GOTO 200
      ENDIF
  200 BESSE1=HELP
      RETURN
      END
      FUNCTION BESSE2(K,N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(9)
      DATA Z,U/0.0D0,1.0D0/
      IF (N.EQ.0) THEN
      IF ((DABS(X).LE.3.75D0).OR.(.NOT.((K.EQ.4).OR.
     1(X.GT.2.D0)))) THEN
      T=X*X/14.0625D0
      A(1)=0.45813D-2
      A(2)=0.360768D-1
      A(3)=0.2659732D0
      A(4)=1.2067492D0
      A(5)=3.0899424D0
      A(6)=3.5156229D0
      A(7)=U
      CALL BESSE3(7,A,T,HELP)
      IF (K.EQ.3) GOTO 200
      T=0.25D0*X*X
      A(1)=0.74D-5
      A(2)=0.1075D-3
      A(3)=0.262698D-2
      A(4)=0.348859D-1
      A(5)=0.23069756D0
      A(6)=0.4227842D0
      A(7)=-0.57721566D0
      CALL BESSE3(7,A,T,R)
      HELP=-HELP*DLOG(0.5D0*X)+R
      GOTO 200
      ENDIF
      IF (K.EQ.3) THEN
      T=3.75D0/X
      A(1)=0.3923767D-2
      A(2)=-0.16476329D-1
      A(3)=0.26355372D-1
      A(4)=-0.20577063D-1
      A(5)=0.9162808D-2
      A(6)=-0.1575649D-2
      A(7)=0.2253187D-2
      A(8)=0.13285917D-1
      A(9)=0.39894228D0
      CALL BESSE3(9,A,T,R)
      HELP=R*DEXP(X)/DSQRT(X)
      GOTO 200
      ENDIF
      T=2.D0/X
      A(1)=0.53208D-3
      A(2)=-0.25154D-2
      A(3)=0.587872D-2
      A(4)=-0.1062446D-1
      A(5)=0.2189568D-1
      A(6)=-0.7832358D-1
      A(7)=1.25331414D0
      CALL BESSE3(7,A,T,R)
      HELP=R/(DSQRT(X)*DEXP(X))
      GOTO 200
      ENDIF
      IF (N.EQ.1) THEN
      IF ((DABS(X).LE.3.75D0).OR.(.NOT.((K.EQ.4).OR.
     1(X.GT.2.D0)))) THEN
      T=X*X/14.0625D0
      A(1)=0.32411D-3
      A(2)=0.301532D-2
      A(3)=0.2658733D-1
      A(4)=0.15084934D0
      A(5)=0.51498869D0
      A(6)=0.87890594D0
      A(7)=0.5D0
      CALL BESSE3(7,A,T,R)
      HELP=X*R
      IF (K.EQ.3) GOTO 200
      T=0.25D0*X*X
      A(1)=-0.4686D-4
      A(2)=-0.110404D-2
      A(3)=-0.1919402D-1
      A(4)=-0.18156897D0
      A(5)=-0.67278579D0
      A(6)=0.15443144D0
      A(7)=U
      CALL BESSE3(7,A,T,R)
      HELP=HELP*DLOG(0.5D0*X)+R/X
      GOTO 200
      ENDIF
      IF (K.EQ.3) THEN
      T=3.75D0/X
      A(1)=-0.4200587D-2
      A(2)=0.17876535D-1
      A(3)=-0.28953121D-1
      A(4)=0.22829673D-1
      A(5)=-0.1031555D-1
      A(6)=0.1638014D-2
      A(7)=-0.3620183D-2
      A(8)=-0.39880242D-1
      A(9)=0.39894228D0
      CALL BESSE3(9,A,T,R)
      HELP=R*DEXP(X)/DSQRT(X)
      GOTO 200
      ENDIF
      T=2.0D0/X
      A(1)=-0.68245D-3
      A(2)=0.325614D-2
      A(3)=-0.780353D-2
      A(4)=0.1504268D-1
      A(5)=-0.365562D-1
      A(6)=0.23498619D0
      A(7)=1.25331414D0
      CALL BESSE3(7,A,T,R)
      HELP=R/(DSQRT(X)*DEXP(X))
      ENDIF
  200 BESSE2=HELP
      RETURN
      END
      SUBROUTINE BESSE3(N,A,T,R)
      DOUBLE PRECISION A,T,R
      DIMENSION A(N)
      DO 10 I=2,N
   10 A(I)=A(I-1)*T+A(I)
      R=A(N)
      RETURN
      END
