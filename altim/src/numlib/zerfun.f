C    ****************************************************************** 
C                                                      66. ZERFUN
C                                                                       
C    PURPOSE                                                            
C    TO CALCULATE ONE ROOT OF AN ARBITRARY COMPLEX FUNCTION F(Z)
C                                                                       
C    USAGE                                                              
C    CALL ZERFUN(Z,H,EP,FUN,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    Z    I Z(2) CONTAINS THE STARTING VALUE Z=Z(1)+IZ(2)
C         O THE RESULTANT ROOT OF F(Z)
C    H    I STEP-LENGTH, SEE REMARK 1
C    EP   I TOLERANCE, SEE REMARK 2
C    FUN  I THE COMPLEX FUNCTION F(Z), SEE REMARK 3
C    IF   O INTEGER EXIT CODE
C           IF=0 NORMAL EXIT
C           IF=1 ABNORMAL EXIT, SEE REMARK 4
C
C    REMARKS
C    1. THE OPTIMAL CHOICE FOR H IS THE APPROXIMATION OF THE
C       DISTANCE BETWEEN THE STARTING VALUE AND THE ROOT.
C    2. THE PROCESS IS STOPPED, WHEN DABS(DZ) < EP HOLDS FOR
C       TWO SUCCESSIVE APPROXIMATIONS.
C    3. THE SUBROUTINE FUN FOR OBTAINING F(Z) MUST BE OF THE
C       FOLLOWING FORM:
C                  SUBROUTINE FUN(Z,F,S)
C                  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                  DIMENSION Z(2),F(2)
C                  <STATEMENTS TO CALCULATE F(1)+IF(2)
C                   FROM Z(1)+IZ(2) AND S=F(1)**2+F(2)**2>
C    4. IF AFTER 200 STEPS NO CONVERGENCE IS OBTAINED FOR THE
C       CALCULATION OF A ROOT, THE PROCESS IS STOPPED
C       WITH IF=1. WHERE Z THEN HAS THE CURRENT VALUE AFTER THESE
C       200 STEPS, THE PROCESS CAN BE RESUMED BY SIMPLY EXECUTING
C       THE PROCEDURE AGAIN.
C
C    METHOD
C    REFER TO RC-TWA-77011/12 OF THE COMPUTER CENTRE
C    OF THE DELFT UNIVERSITY OF TECHNOLOGY
C
C    *******************************************************************
      SUBROUTINE ZERFUN(Z,H,EP,FUN,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FUN
      INTEGER W
      LOGICAL BOL
      DIMENSION Z(2),F(2),Y(2),V(2),U(2)
      DATA ZZ,UU/0.D0,1.D0/
      M=0
      W=0
      IF=0
      D=2.D-3
      V(1)=Z(1)
      V(2)=Z(2)
      Y(1)=ZZ
      Y(2)=ZZ
   10 M=M+1
      IF (M.GE.200) THEN
      IF=1
      GOTO 100
      ENDIF
      BOL=.FALSE.
      CALL FUN(Z,F,F1)
      S1=F(1)
      T1=F(2)
      DO 20 I=1,2
      R=Z(I)-V(I)
      IF (Y(I)*R.LT.ZZ) W=W+1
      Y(I)=R
   20 V(I)=Z(I)
      IF (W.GT.3) THEN
      W=0
      H=H*0.5D0
      ENDIF
   30 IF (BOL) THEN
      S=S-S1
      T=T-T1
      R=S*S+T*T
      S2=(S1*S+T1*T)/R
      T2=(T1*S-S1*T)/R
      R=V1*S2-V2*T2
      V2=V1*T2+V2*S2
      V1=R
      S=S1
      T=T1
      Z(1)=Z(1)+V1
      Z(2)=Z(2)+V2
      R=DSQRT(V1*V1+V2*V2)
      IF (R.LT.EP) GOTO 100
      CALL FUN(Z,F,F0)
      IF (F0.GE.F1) THEN
      Z(1)=Z(1)-V1
      Z(2)=Z(2)-V2
      GOTO 40
      ENDIF
      S1=F(1)
      T1=F(2)
      F1=F0
      GOTO 30
      ENDIF
   40 N=0
      K=8
      BOL=.FALSE.
      S=S1
      T=T1
      F0=F1
      IF (F0.LT.1.D-50) GOTO 100
      U(2)=Z(2)+D
      U(1)=Z(1)
      CALL FUN(U,F,F2)
      S2=F(1)
      T2=F(2)
      U(1)=Z(1)+D
      U(2)=Z(2)
      CALL FUN(U,F,F1)
      S1=F(1)
      T1=F(2)
      V1=S1+S2-S-S
      V2=T1+T2-T-T
      S2=S1-S2
      T2=T1-T2
      S1=(V2+V1)/D
      T1=(V2-V1)/D
      S2=(S2-V2)/(D*D)
      T2=(T2+V1)/(D*D)
      F1=0.5D0*(S1*S1-T1*T1)
      F2=S1*T1
      R=F1-4.D0*(S*S2-T*T2)
      T2=F2-4.D0*(T*S2+S*T2)
      S2=R
      R1=S*S1+T*T1
      R2=T*S1-S*T1
      R=S1*S1+T1*T1
      IF (R.LT.1.D-50) GOTO 50
      R=2.D0/R
      V1=R1*R
      V2=R2*R
      Q=V1*V1+V2*V2
      IF (D.GT.Q) D=Q
      Q=DSQRT(Q)
      IF (Q.GT.5.D-3) GOTO 50
      C1=S2+F1
      C2=T2+F2
      R=C1*C1+C2*C2
      IF (R.LT.1.D-50) THEN
      V1=-V1
      V2=-V2
      GOTO 90
      ENDIF
      R=2.D0/R
      Q=(F1*C1+F2*C2)*R
      F2=(F2*C1-F1*C2)*R
      F1=Q
      Q=DSQRT(F1*F1+F2*F2)
      IF (Q.LT.0.95D0) GOTO 50
      IF (Q.GT.1.05D0) THEN
      R=V1*F2+V2*F1
      V1=V2*F2-V1*F1
      V2=-R
      GOTO 60
      ENDIF
      BOL=.TRUE.
   50 R=DSQRT(S2*S2+T2*T2)
      V1=S2
      S2=DSQRT(R+DABS(S2))
      T2=T2/S2
      IF (V1.LT.ZZ) THEN
      Q=S2
      S2=T2
      T2=Q
      ENDIF
      Q=DSIGN(UU,S1*S2+T1*T2)
      IF (Q.EQ.ZZ) Q=UU
      S1=S1+Q*S2
      T1=T1+Q*T2
      Q=-4.D0/(S1*S1+T1*T1)
      V1=Q*(S*S1+T*T1)
      V2=Q*(T*S1-S*T1)
   60 Q=DABS(V1)+DABS(V2)
      IF (Q.LT.EP) GOTO 90
      U(1)=Z(1)+V1
      U(2)=Z(2)+V2
      CALL FUN(U,F,F1)
      S1=F(1)
      T1=F(2)
      IF (F1.GT.F0) GOTO 80
   70 N=N+1
      IF ((N.EQ.10).AND.(K.GT.0)) THEN
      Q=H/DSQRT(R1*R1+R2*R2)
      V1=-R1*Q
      V2=0.5D0*(V1+R2*Q)
      Z(1)=Z(1)+V2-V1
      Z(2)=Z(2)+V2
      GOTO 10
      ENDIF
      S2=V1+V1
      T2=V2+V2
      U(1)=Z(1)+S2
      U(2)=Z(2)+T2
      CALL FUN(U,F,F2)
      IF (F2.GT.F1) THEN
      Z(1)=Z(1)+V1
      Z(2)=Z(2)+V2
      GOTO 30
      ENDIF
      S1=F(1)
      T1=F(2)
      V1=S2
      V2=T2
      F1=F2
      GOTO 70
   80 F2=F1
      Q=0.5D0*Q
      S2=S1
      T2=T1
      K=K-1
      IF (Q.LT.EP) THEN
      IF (K.LT.0) GOTO 100
      K=-1
      V1=-4.D0*V1
      V2=-4.D0*V2
      GOTO 60
      ENDIF
      V1=0.5D0*V1
      V2=0.5D0*V2
      U(1)=Z(1)+V1
      U(2)=Z(2)+V2
      CALL FUN(U,F,F1)
      S1=F(1)
      T1=F(2)
      IF (F1.GE.F0) GOTO 80
      Z(1)=U(1)
      Z(2)=U(2)
      GOTO 30
   90 Z(1)=Z(1)+V1
      Z(2)=Z(2)+V2
  100 CONTINUE
      RETURN
      END
