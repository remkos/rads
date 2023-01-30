C    ****************************************************************** 
C                                                      13. FUNAFG
C                                                                       
C    PURPOSE                                                            
C    TO CALCULATE THE FIRST OR THE SECOND DERIVATIVE
C    OF A REAL FUNCTION, USING ROMBERG ITERATION
C                                                                       
C    USAGE                                                              
C    CALL FUNAFG(A,X,H,HM,N,M,FUN,P,G,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    O THE RESULTANT VALUE
C    X    I THE ABCISSA OF THE DERIVATIVE
C    H    I INITIAL STEP, SEE REMARK 1
C         O THE OBTAINED STARTING STEP H, SEE REMARKS
C    HM   I MAXIMAL STEP H
C         O APPROXIMATION OF THE ABSOLUTE ERRROR IN A
C    N    I N=1  A=FIRST  DERIVATIVE
C         I N=2  A=SECOND DERIVATIVE
C    M    I MAXIMAL NUMBER OF ROMBERG STEPS
C    FUN  I GIVEN FUNCTION
C    P      P(M) IS AN AUXILIARY INTEGER WORKING VECTOR
C    G      G(M) IS AN AUXILIARY INTEGER WORKING VECTOR
C    IF   I    IF=0 FROM THE STARTING H IT IS TRIED
C                   TO FIND A OPTIMAL ONE FIRST
C              IF=1 THE GIVEN VALUE OF H IS USED
C         O INTEGER ERROR CODE
C              IF=0 NORMAL EXIT
C              IF=1 NUMBER OF STEPS=M, WITHOUT FINDING
C                   AN OPTIMAL VALUE FOR A
C
C    SUBROUTINES REQUIRED
C    1 FUNCTION FUN(X)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      <STATEMENTS YIELDING THE FUNCTION VALUE F(X)
C    2 FUNAFG USES SUBROUTINE FUNAF1
C
C    REMARKS
C    1. THE STEP H MUST ALWAYS BE SMALLER THAN THE RADIUS OF
C       CONVERGENCE OF THE TAYLOR SERIES. THIS NEED NOT BE
C       KNOWN EXACTLY. A REASONABLE BIG VALUE WILL SATISFY
C    2. IF IF=1 IN THE OUTPUT, A CAN STILL BE A FAIRLY GOOD
C       APPROXIMATION, FOR WHICH HM CAN SERVE AS A MEASURE.
C       ONE CAN ALWAYS START THE PROGRAM AGAIN, WITH A GREATER
C       VALUE FOR M AND IF=1 WITH THE RESULTANT VALUE FOR H.
C
C    METHOD
C    REFER TO RC-TWA-81001 OF COMPUTER CENTRE D.U.T.
C
C    *******************************************************************
      SUBROUTINE FUNAFG(A,X,H,HM,N,M,FUN,P,G,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(M),G(M)
      EXTERNAL FUN
      W=0.5D0
      L=1
      K=1
      S=H
      HS=H
      FX=FUN(X)
      F=W*DABS(FX)
      CALL FUNAF1(N,S,T,U,V,X,F,FX,FUN)
      G(1)=U
      IF (IF.EQ.1) GOTO 30
      T=DABS(T)
      IF (T.LT.V*0.2D0) THEN
      V1=-1.D0
      V2=0.2D0
      W=1.D0/W
      ELSE
      IF (T.LT.V*2.5D0) THEN
      GOTO 30
      ELSE
      V1=1.D0
      V2=2.D0
      K=0
      ENDIF
      ENDIF
      V2=V2*V1
   10 L=L+K
      IF (L.GT.M) THEN
      L=M
      HS=HS*W
      DO 20 I=2,L
   20 G(I-1)=G(I)
      ENDIF
      S=S*W
      CALL FUNAF1(N,S,T,U,V,X,F,FX,FUN)
      G(L)=U
      IF ((V1*DABS(T).GT.V2*V).AND.(S*W.LT.HM)) GOTO 10
   30 V1=S
      IF=0
      IF (W.GT.1.D0) THEN
      W=1.D0/W
      S=HS
      J=(L)/2
      DO 40 I=1,J
      T=G(I)
      G(I)=G(L-I+1)
   40 G(L-I+1)=T
      ENDIF
      P(1)=1.D0
      R=W*W
      TM=1.D70
      F=F*1.D-15
      DO 60 I=2,M
      P(I)=P(I-1)*R
      IF (I.GT.L) THEN
      S=S*W
      L=I
      CALL FUNAF1(N,S,T,U,V,X,F,FX,FUN)
      G(I)=U
      ENDIF
      A=G(1)
      DO 50 J=I-1,1,-1
      T=P(I-J+1)
   50 G(J)=(T*G(J)-G(J+1))/(T-1.D0)
      T=DABS(G(1)-A)
      IF (((T.GE.TM).OR.(T.LT.F)).AND.(I.GT.2)) GOTO 70
   60 TM=T
      IF=1
   70 HM=TM
      H=V1
      RETURN
      END
      SUBROUTINE FUNAF1(N,S,T,U,V,X,F,FX,FUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FUN
      Y=FUN(X+S)
      Z=FUN(X-S)
      IF (N.EQ.1) THEN
      T=Y-Z
      U=0.5D0*T/S
      V=0.5D0*(DABS(Y)+DABS(Z))
      ELSE
      T=Y+Z-2.D0*FX
      U=T/(S*S)
      V=0.25D0*(DABS(Y)+DABS(Z))+F
      ENDIF
      RETURN
      END
