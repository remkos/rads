C    ****************************************************************** 
C                                                      56. POLBAR
C                                                                       
C    PURPOSE                                                            
C    TO CALCULATE A CERTAIN NUMBER OF ROOTS OF A POLYNOMIAL
C    WITH REAL COEFFICIENTS BY BAIRSTOW'S METHOD.
C                                                                       
C    USAGE                                                              
C    CALL POLBAR(N1,M,IE,A,E,B,P,Q,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I DEGREE OF THE POLYNOMIAL
C         O DEGREE OF DEFLATED POLYNOMIAL
C    M    I NUMBER OF WANTED ROOTS
C         O NUMBER OF CALCULATED ROOTS, SEE REMARK 5
C    IE   I FIRST DIMENSION OF E AS DECLARED IN THE
C           DIMENSION STATEMENT OF THE CALLING PROGRAM
C    A    I ARRAY A(N1) CONTAINS THE COEFFICIENTS
C           A(I-1) FOR I=1(1)N1, WHERE N1=N+1
C         O COEFFICIENTS OF DEFLATED POLYNOMIAL
C    E    O ARRAY E(N,2) CONTAINS THE CALCULATED ROOTS, SEE REMARK 1
C    B      AUXILIARY WORKING VECTOR OF AT LEAST 4*N+19 ELEMENTS.
C    P,Q  I STARTING VALUES OF THE COEFFICIENTS
C           OF THE QUADRATIC FACTOR X**2+P*X+Q
C    EP   I TOLERANCE, SEE REMARK 2
C    IF   O INTEGER EXIT CODE
C           IF=0 NORMAL EXIT
C           IF=1 ABNORMAL EXIT, SEE REMARK 3
C
C    REMARKS
C    1. THE CALCULATED ROOTS WILL APPEAR IN REVERSED ORDER IN
C       THE MATRIX E. SO SO THE FIRST ROOT WILL BE E(N,1)+I*E(N,2).
C       IF M ROOTS ARE CALCULATED, THE VALUES E(J,1)+I*E(J,2)
C       MUST BE PUT INTO OUTPUT FOR J=N-M+1/N
C    2. THE PROCESS FOR THE CALCULATION OF THE QUADRATIC FACTOR
C       WILL BE STOPPED IF DABS((DP)**2+(DQ)**2) < EP.
C    3. IF AFTER 200 STEPS NO CONVERGENCE IS OBTAINED FOR THE
C       CALCULATION OF A QUADRATIC FACTOR, THE PROCESS IS STOPPED
C       WITH IF=1. IT IS THEN POSSIBLE TO TRY TO OBTAIN THE
C       REMAINING ROOTS FROM THE DEFLATED POLYNOMIAL A OF DEGREE
C       N WITH DIFFERENT VALUES FOR P AND Q.
C    4. IF ALREADY APPROXIMATIONS EXIST OF THE ROOTS AND THUS OF
C       P AND Q, THEN THE SUBROUTINE CAN BE USED AS FOLLOWS
C             K=N/2
C             DO 20 J=1,K
C          20 CALL POLBAR(N,2,IE,A,E,B,P(J),Q(J),EP,IF)
C    5. IF M IS GIVEN AN ODD VALUE, THEN M IS SET TO M+1.
C    6. IF M>=N-2 THEN M IS SET TO N
C    7. A(N1) MUST HAVE THE VALUE 1.
C    8. POLBAR USES SUBROUTINES POLBA1 AND POLBA2
C
C    METHOD
C    A GENERALISED BAIRSTOW METHOD, REFER TO RC-TWA-77011/12
C    OF THE COMPUTER CENTRE OF THE DELFT UNIVERSITY OF TECHNOLOGY
C
C    *******************************************************************
      SUBROUTINE POLBAR(N,M,IE,A,E,B,P,Q,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BOL,BOL1
      DIMENSION C(4),A(*),E(IE,2),B(*)
      DATA ZZ,UU/0.D0,1.D0/
      K=N
      IF=0
      KB=N+6
      KC=KB+KB
      KD=KC+KB
      IF (M.GT.N) M=N
      DO 4 I=1,N+1
    4 B(I)=A(I)
   10 I=N-K
      IF (I.GE.M) GOTO 160
      R1=UU
      R2=UU
      J=0
      X=P
      Y=Q
      L=0
      DO 20 I=-4,K+1
      B(I+KB)=ZZ
      B(I+KC)=ZZ
   20 B(I+KD)=ZZ
      B(K+KB)=UU
      B(K+KC-1)=UU
      B(K+KD-2)=UU
      BOL=.TRUE.
      BOL1=.FALSE.
   30 IF (BOL1) THEN
      L=L+1
      X=P
      Y=Q
      IF (L.EQ.K-1) L=0
      ENDIF
      CALL POLBA1(X,Y,B,B0,B1,F,K-1,L,KB,1)
   40 IF (F.LT.1.D-50) GOTO 130
      IF (DSQRT(R1*R1+R2*R2).LT.EP) GOTO 130
      J=J+1
      IF (J.GT.200) THEN
      IF=1
      GOTO 160
      ENDIF
      DO 50 I=L,K-2
   50 B(I+KB)=B(I+KB+1)
      CALL POLBA1(X,Y,B,R,S,T,K-2,L,KC,KB)
      DO 60 I=0,2
   60 C(I+1)=B(I+L+KC-3)-B(I+L+KC)
      C1=C(2)
      C2=C(3)
      C0=C(1)+B1+X*C1
      C3=C1+X*C2
      IF (BOL) THEN
      R1=C2*B0-C3*B1
      R2=C0*B1-C1*B0
      R=C0*C2-C1*C3
      S=X
      T=Y
      IF (DABS(R).LT.1.D-10) GOTO 70
      R1=-R1/R
      R2=-R2/R
      IF (R1*R1+R2*R2.GT.UU) GOTO 70
      X=X+R1
      Y=Y+R2
      CALL POLBA1(X,Y,B,A11,A22,F1,K-1,L,KB,1)
      IF (F1.GE.F) GOTO 70
      B0=A11
      B1=A22
      F=F1
      GOTO 40
   70 X=S
      Y=T
      BOL=.FALSE.
      ENDIF
      DO 80 I=L,K-2
   80 B(I+KC)=B(I+KC+1)
      CALL POLBA1(X,Y,B,R,S,T,K-3,L,KD,KC)
      DO 90 I=0,3
   90 C(I+1)=(B(I+L+KD)-B(I+L+KD-4))*2.D0
      A11=C0*C0+C1*C1
      A22=C3*C3+C2*C2
      A12=C0*C3+C1*C2
      R=B1+X*B0
      A11=A11+B0*(C(1)+2.D0*C1)+R*C(2)
      A12=A12+B0*(C(2)+C2)+R*C(3)
      A22=A22+B0*C(3)+R*C(4)
      S=A11-A22
      T=A11+A22
      R=DSQRT(S*S+4.D0*A12*A12)
      S=0.5D0*(T+R)
      T=S-R
      F1=C0*B0+C1*B1
      F2=C3*B0+C2*B1
      IF (T.LT.ZZ) THEN
      IF (A11.GT.A22) THEN
      R1=-A12
      R2=A11-T
      ELSE
      R1=A22-T
      R2=-A12
      ENDIF
      R=(R1*F1+R2*F2)/F
      T=2.D0*T/F
      R=DSIGN(UU,R)*2.D0/(DSQRT(R*R-T*(R1*R1+R2*R2))+DABS(R))
      R1=-R*R1
      R2=-R*R2
      BOL1=.FALSE.
      GOTO 100
      ENDIF
      BOL=(F.LT.5.D-4*T)
      BOL1=(F1*F1+F2*F2.LT.F*1.D-3)
      R1=A22*F1-A12*F2
      R2=A11*F2-A12*F1
      R=A11*A22-A12*A12
      IF (DABS(R).LT.1.D-10) THEN
      R=F1*F1+F2*F2
      IF (R.LT.1.D-30) THEN
      BOL=.FALSE.
      GOTO 30
      ENDIF
      R=F/R
      R1=-R*F1
      R2=-R*F2
      GOTO 100
      ENDIF
      R1=-R1/R
      R2=-R2/R
  100 S=X
      T=Y
      X=S+R1
      Y=T+R2
      CALL POLBA1(X,Y,B,B0,B1,F1,K-1,L,KB,1)
      R=DSQRT(R1*R1+R2*R2)
      IF (F1.GT.F) GOTO 120
  110 X=X+R1
      Y=Y+R2
      CALL POLBA1(X,Y,B,B0,B1,F2,K-1,L,KB,1)
      IF (F2.LT.F1) THEN
      R1=R1+R1
      R2=R2+R2
      F1=F2
      GOTO 110
      ENDIF
      X=S+R1
      Y=T+R2
      GOTO 30
  120 R1=0.5D0*R1
      R2=0.5D0*R2
      R=0.5D0*R
      F2=F1
      IF (R.LT.EP) GOTO 130
      X=S+R1
      Y=T+R2
      CALL POLBA1(X,Y,B,B0,B1,F1,K-1,L,KB,1)
      IF (F1.GT.F) GOTO 120
      GOTO 30
  130 B(K)=ZZ
      B(K+1)=ZZ
      CALL POLBA2(X,Y,K,E,IE)
      K=K-2
      IF (K.GT.2) THEN
      DO 140 I=L-1,0,-1
  140 B(I+1)=B(I+KB)
      DO 150 I=L,K
  150 B(I+1)=B(I+KB+2)
      GOTO 10
      ENDIF
      IF (K.EQ.2) THEN
      X=B(3+KB)
      Y=B(2+KB)
      CALL POLBA2(X,Y,K,E,IE)
      ELSE
      E(1,2)=ZZ
      E(1,1)=-B(2+KB)
      ENDIF
      K=0
      B(1)=ZZ
      B(2)=ZZ
      B(3)=ZZ
  160 M=N-K
      DO 170 I=1,N+1
  170 A(I)=B(I)
      N=K
      RETURN
      END
      SUBROUTINE POLBA1(X,Y,B,R,S,T,K,L,KW,KZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(*)
      S=1.D0/Y
      DO 10 I=0,L-1
   10 B(I+KW)=S*(B(I+KZ)-B(I+KW-2)-X*B(I+KW-1))
      DO 20 I=K,L,-1
   20 B(I+KW)=B(I+KZ)-X*B(I+KW+1)-Y*B(I+KW+2)
      S=B(L+KW+1)-B(L+KW-1)
      R=B(L+KW)-B(L+KW-2)+X*S
      T=0.5D0*(S*S+R*R)
      RETURN
      END
      SUBROUTINE POLBA2(X,Y,K,E,IE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(IE,2)
      R=X*X-4.D0*Y
      J=K-1
      S=0.5D0*DSQRT(DABS(R))
      IF (R.LT.0.D0) THEN
      E(K,1)=-0.5D0*X
      E(J,1)=E(K,1)
      E(K,2)=-S
      E(J,2)=S
      ELSE
      E(K,2)=0.D0
      E(J,2)=0.D0
      T=DABS(X)
      IF (T+DABS(Y).LT.1.D-40) THEN
      E(K,1)=0.D0
      E(J,1)=0.D0
      ELSE
      R=-DSIGN(1.D0,X)*(0.5D0*T+S)
      E(K,1)=R
      E(J,1)=Y/R
      ENDIF
      ENDIF
      RETURN
      END
