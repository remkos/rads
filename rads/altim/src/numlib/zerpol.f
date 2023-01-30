C    ****************************************************************** 
C                                                      67. ZERPOL
C                                                                       
C    PURPOSE                                                            
C    TO CALCULATE A CERTAIN NUMBER OF ROOTS OF A POLYNOMIAL
C    PN(X) WITH REAL OR COMPLEX COEFFICIENTS.
C                                                                       
C    USAGE                                                              
C    CALL ZERPOL(N1,M,IA,IE,A,E,B,P,Q,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N1   I DEGREE (PLUS 1) OF THE POLYNOMIAL
C           SUM(I=1:N1) A(I)*X**(I-1)
C           NOTE N=N1-1 IS THE DEGREE OF PN(X)=0
C         O DEGREE (PLUS 1) OF DEFLATED POLYNOMIAL
C    M    I NUMBER OF WANTED ROOTS, M<=N
C         O NUMBER OF CALCULATED ROOTS
C    IA   I FIRST DIMENSION OF A AND B AS DECLARED IN THE
C           DIMENSION STATEMENT OF THE CALLING PROGRAM
C    IE   I FIRST DIMENSION OF E AS DECLARED IN THE
C           DIMENSION STATEMENT OF THE CALLING PROGRAM
C    A    I ARRAY A(N1,2) CONTAINS THE COEFFICIENTS
C           A(I,1)+IA(I,2) OF PN(X) FOR I=1:N1
C         O COEFFICIENTS OF DEFLATED POLYNOMIAL
C    E    O ARRAY E(N,2) CONTAINS THE CALCULATED ROOTS, SEE REMARK 1
C    B      AUXILIARY WORKING ARRAY OF DIMENSION AT LEAST (N1,8)
C    P,Q  I STARTING VALUE FOR Z0 IS P+IQ
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
C    2. THE PROCESS FOR THE CALCULATION OF A ROOT IS STOPPED, WHEN
C       DABS(DZ) < EP HOLDS FOR TWO SUCCESSIVE APPROXIMATIONS.
C    3. IF AFTER 100 STEPS NO CONVERGENCE IS OBTAINED FOR THE
C       CALCULATION OF A ROOT, THE PROCESS IS STOPPED WITH IF=1.
C       IT IS THEN ALWAYS POSSIBLE TO TRY TO OBTAIN THE
C       REMAINING ROOTS FROM THE DEFLATED POLYNOMIAL A OF DEGREE
C       N WITH DIFFERENT VALUES FOR P AND Q.
C    4. IF ALREADY APPROXIMATIONS EXIST OF THE ROOTS AND THUS OF
C       P AND Q, THEN THE SUBROUTINE CAN BE USED AS FOLLOWS
C             DO 20 J=1,M
C          20 CALL ZERPOL(N,1,IA,IE,A,E,B,P(J),Q(J),EP,IF)
C    5. A(N1) MUST HAVE THE VALUE 1.
C    6. ZERPOL USES THE SUBROUTINES ZERPO1, ZERPO2 AND ZERPO3.
C
C    METHOD
C    A GENERALISED LAQUERRE METHOD, REFER TO RC-TWA-77011/12
C    OF THE COMPUTER CENTRE OF THE DELFT UNIVERSITY OF TECHNOLOGY
C
C    *******************************************************************
      SUBROUTINE ZERPOL(N1,M,IA,IE,A,E,B,P,Q,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BOL
      DIMENSION A(IA,2),E(IE,2),B(IA,8)
      DATA ZZ,UU/0.D0,1.D0/
      N=N1-1
      X=A(N1,1)
      Y=A(N1,2)
      T=UU/(X*X+Y*Y)
      B(N1,1)=UU
      B(N1,2)=ZZ
      DO 10 K=N,1,-1
      B(K,1)=(A(K,1)*X+A(K,2)*Y)*T
   10 B(K,2)=(A(K,2)*X-A(K,1)*Y)*T
      IF (M.EQ.0) GOTO 80
      IF (M.GT.N) M=N
      T=EP*EP
      X=B(1,1)**2+B(1,2)**2
      DO 20 K=2,N1
      Y=B(K,1)**2+B(K,2)**2
      IF (X.GT.Y*T) THEN
      L=K
      GOTO 30
      ENDIF
   20 X=Y
   30 DO 40 K=N-L+3,N
      E(K,1)=ZZ
   40 E(K,2)=ZZ
      L=L-2
      IF (L.GT.0) THEN
      DO 50 K=L+1,N1
      B(K-L,1)=B(K,1)
   50 B(K-L,2)=B(K,2)
      ENDIF
      IF=0
      L=N1-L
      Y=ZZ
      BOL=.TRUE.
      DO 60 K=L,2,-1
      IF (K.EQ.N1-M) GOTO 80
      IF ((DABS(Y).GT.EP).AND.(BOL)) THEN
      Y=-Y
      BOL=.FALSE.
      ELSE
      X=P
      Y=Q
      BOL=.TRUE.
      ENDIF
      IF (K.EQ.2) THEN
      X=-B(1,1)
      Y=-B(1,2)
      ELSE
      CALL ZERPO1(K,IA,B,B,X,Y,EP,IF)
      IF (IF.EQ.1) GOTO 70
      ENDIF
      S=X
      T=Y
      CALL ZERPO1(N1,IA,A,B,X,Y,EP,IF)
      IF (IF.EQ.1) GOTO 70
      E(K-1,1)=X
      E(K-1,2)=Y
      CALL ZERPO3(K,IA,B,S,T)
   60 CONTINUE
      GOTO 80
   70 M=N1-K
   80 CONTINUE
      DO 90 K=N1-M,1,-1
      A(K,1)=B(K,1)
   90 A(K,2)=B(K,2)
      DO 100 K=N1-M+1,N1
      A(K,1)=ZZ
  100 A(K,2)=ZZ
      N1=N1-M
      END
      SUBROUTINE ZERPO1(N,IA,A,B,X,Y,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BOL1,BOL2,BOL3,BOL4
      DIMENSION A(IA,2),B(IA,8)
      DATA ZZ,UU/0.D0,1.D0/
      BOL1=.FALSE.
      BOL2=.FALSE.
      BOL3=.FALSE.
      BOL4=.FALSE.
      DO 10 J=1,2
      DO 10 I=1,N
      T=A(I,J)
      B(I,J*3)=T
      T=T*(I-1)
      B(I,J*3+1)=T
   10 B(I,J*3+2)=T*(I-2)
      M=0
      C=10**(35./N)
      CALL ZERPO2(N,IA,0,B,X,Y,S,T)
      F=S*S+T*T
      G=A(1,1)**2+A(1,2)**2
      L=1
   20 IF (F.LT.G*1.D-29) GOTO 110
      M=M+1
      IF (M.EQ.100) THEN
      IF=1
      GOTO 110
      ENDIF
      IF (BOL1) THEN
      IF (BOL2) GOTO 30
      R1=R1-S
      R2=R2-T
      BOL2=.TRUE.
   30 R1=R1-S
      R2=R2-T
      D=R1*R1+R2*R2
      IF (D.EQ.ZZ) GOTO 110
      F2=(S*R1+T*R2)/D
      R2=(T*R1-S*R2)/D
      R1=V1*F2-V2*R2
      V2=V2*F2+V1*R2
      V1=R1
      R1=S
      R2=T
      GOTO 70
      ENDIF
      R1=S
      R2=T
      CALL ZERPO2(N,IA,1,B,X,Y,S1,T1)
      F1=4.D0/(DABS(S)+DABS(T)+DABS(S1)+DABS(T1))
      S=S*F1
      T=T*F1
      S1=S1*F1
      T1=T1*F1
      D=S1*S1+T1*T1
      IF (D.EQ.ZZ) GOTO 110
      V1=-(S*S1+T*T1)/D
      V2=-(T*S1-S*T1)/D
      F2=5.D-3
      D=DSQRT(V1*V1+V2*V2)
      BOL3=(D*L.LT.F2)
      IF (BOL3.AND.BOL4) GOTO 60
   40 F2=N
      BOL4=BOL3
      CALL ZERPO2(N,IA,2,B,X,Y,S2,T2)
      S2=S2*F1
      T2=T2*F1
      D=S*S2-T*T2
      T2=S*T2+T*S2
      D2=D
      V1=S1*S1-T1*T1
      V2=2.D0*S1*T1
      F1=V1-S2
      U=V2-T2
      D=F1*F1+U*U
      IF (D.EQ.ZZ) THEN
      U=UU
      GOTO 50
      ENDIF
      U=DSQRT((V1*V1+V2*V2)/D)
   50 L=U+0.5D0
      F1=F2/(F2-UU)
      V1=0.5D0*(V1-F1*S2)
      V2=0.5D0*(V2-F1*T2)
      S2=DSQRT(V1*V1+V2*V2)
      S2=DSQRT(S2+DABS(V1))
      T2=V2/S2
      IF (V1.LT.ZZ) THEN
      D=S2
      S2=T2
      T2=D
      ENDIF
      D=DSIGN(UU,S1*S2+T1*T2)
      IF (D.EQ.ZZ) D=UU
      S1=S1/(F2-UU)+D*S2
      T1=T1/(F2-UU)+D*T2
      D=-F1/(S1*S1+T1*T1)
      V1=D*(S*S1+T*T1)
      V2=D*(T*S1-S*T1)
      BOL3=.FALSE.
      GOTO 70
   60 BOL1=(L.EQ.1)
      IF (L.LT.1) L=1
      V1=L*V1
      V2=L*V2
   70 S2=X
      T2=Y
   80 X=S2+V1
      Y=T2+V2
      IF (V1*V1+V2*V2.GT.C) THEN
      V1=0.1D0*V1
      V2=0.1D0*V2
      GOTO 80
      ENDIF
      IF (DABS(V1)+DABS(V2).LT.EP) GOTO 110
      CALL ZERPO2(N,IA,0,B,X,Y,S,T)
      F1=S*S+T*T
      IF (F1.LT.G*1.D-29) GOTO 110
   90 IF (F1.LT.F) GOTO 100
      V1=0.1D0*V1
      V2=0.1D0*V2
      IF (DABS(V1)+DABS(V2).LT.EP) THEN
      X=S2
      Y=T2
      GOTO 110
      ENDIF
      X=S2+V1
      Y=T2+V2
      F2=F1
      CALL ZERPO2(N,IA,0,B,X,Y,S,T)
      F1=S*S+T*T
      GOTO 90
  100 X=X+V1
      Y=Y+V2
      CALL ZERPO2(N,IA,0,B,X,Y,S1,T1)
      F2=S1*S1+T1*T1
      IF (F1.GE.F2) THEN
      V1=V1+V1
      V2=V2+V2
      S=S1
      T=T1
      F1=F2
      GOTO 100
      ENDIF
      X=S2+V1
      Y=T2+V2
      F=F1
      GOTO 20
  110 CONTINUE
      RETURN
      END
      SUBROUTINE ZERPO2(N,IA,K,B,X,Y,S,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(IA,8)
      L=K+3
      J=K+6
      T=B(N,J)
      S=B(N,L)
      DO 10 I=N-1,K+1,-1
      R=X*S-Y*T+B(I,L)
      T=X*T+Y*S+B(I,J)
   10 S=R
      RETURN
      END
      SUBROUTINE ZERPO3(N,IA,B,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(IA,8)
      S=B(N,1)
      T=B(N,2)
      DO 10 J=N-1,2,-1
      R=X*S-Y*T+B(J,1)
      T=X*T+Y*S+B(J,2)
      S=R
      B(J,1)=S
   10 B(J,2)=T
      DO 20 J=2,N
      B(J-1,1)=B(J,1)
   20 B(J-1,2)=B(J,2)
      RETURN
      END
