C    ****************************************************************** 
C                                                      51. MINGRF
C                                                                       
C    PURPOSE                                                            
C    TO CALCULATE A ROOT OF N NON-LINEAR EQUATIONS IN N UNKNOWNS
C    WITH THE GENERALISED REGULA FALSI METHOD
C                                                                       
C    USAGE                                                              
C    CALL MINGRF(N,M,H,X,EP,EP1,A,B,C,D,F,V,W,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I NUMBER OF VARIABLES
C    M    I M=N+1
C    H    I INITIAL STEP, SEE REMARK 1
C    X    I X(N) CONTAINS THE INITIAL VALUES OF THE VARIABLES
C         O X(N) THE FOUND ROOT OR THE LAST APPROXIMATION (IF=2)
C    EP   I EP AND EP1 ARE TOLERANCES, SEE IF AND REMARK 2
C    EP1  I
C    G    I SUBROUTINE YIELDING THE N NON-LINEAR EQUATIONS,
C           SEE UNDER SUBROUTINES REQUIRED
C    A      A(M,M) IS AN AUXILIARY REAL WORKING ARRAY
C    B      B(M,M) IS AN AUXILIARY REAL WORKING ARRAY
C    C      C(M) IS AN AUXILIARY INTEGER WORKING VECTOR
C    D      D(M) IS AN AUXILIARY INTEGER WORKING VECTOR
C    F    O F(N) CONTAINS THE FUNCTION VALUES IN THE RESULT
C    V      V(M) IS AN AUXILIARY REAL WORKING VECTOR
C    W      W(M) IS AN AUXILIARY REAL WORKING VECTOR
C    IF   O INTEGER ERROR CODE
C              IF=0 F(X)**2 < EP1 , A ROOT IS FOUND
C              IF=1 DX**2 < EP    , A ROOT IS FOUND
C              IF=2 NUMBER OF STEPS > 100
C              IF=3 MATRIX SINGULAR, TRY OTHER INITIAL VALUE
C
C    SUBROUTINES REQUIRED
C    1 SUBROUTINE G(N,X,F)
C      DOUBLE PRECISION X,F
C      DIMENSION X(N),F(N)
C      <STATEMENTS YIELDING THE N NON-LINEAR VALUES F(X)>
C    2 MINGRF USES SUBROUTINE INVGAU
C
C    REMARKS
C    1. A MATRIX IS FILLED WITH THE VECTORS X+H*E(J), J=1/N
C       WHERE E(J) IS THE JTH UNIT VECTOR. IT WOULD BE ADVANTEOUS
C       IF THIS SPACE CONTAINS THE ROOT.
C    2. IN THE EXIT CONDITIONS THE SQUARES OF THE LENGTHS OF F
C       AND DX ARE USED, SO EP AND EP1 CAN BE VERY SMALL.
C    3. USE MINGRF FOR SMALL VALUES OF N.
C
C    METHOD
C    REFER TO RC-TWA-74002 OF COMPUTER CENTRE D.U.T.
C
C    *******************************************************************
      SUBROUTINE MINGRF(N,M,H,X,EP,EP1,G,A,B,C,D,F,V,W,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER C,D,Q,Q1
      DIMENSION X(N),A(M,M),B(M,M),C(M),D(M),F(N),V(M),W(M)
      EXTERNAL G
      DATA Z/0.D0/
      IF=0
      P=1.D70
      Q=0
      K1=M
   10 Q1=0
      DO 40 J=1,M
      IF (J.EQ.M) GOTO 20
      S=X(J)
      X(J)=S+H
   20 CALL G(N,X,F)
      T=Z
      DO 30 I=1,N
      A(I,J)=X(I)
      B(I,J)=F(I)
   30 T=T+F(I)*F(I)
      IF (J.LT.M) X(J)=S
      B(M,J)=1.D0
      IF (T.LT.P) THEN
      P=T
      L=J
      ENDIF
   40 CONTINUE
      CALL INVGAU(M,M,B,C,D,V,W,T,EP,IF)
      IF (IF.EQ.1) THEN
      IF=3
      GOTO 180
      ENDIF
      T=Z
   50 Q1=Q1+1
      IF (Q1.EQ.10) THEN
      H=0.1*H
      T=0.1*DSQRT(T)
      IF (H.LT.T) H=T
      IF (H.LT.0.001) GOTO 50
      Q1=0
      GOTO 10
      ENDIF
      Q=Q+1
      IF (Q.GT.100) THEN
      IF=2
      GOTO 160
      ENDIF
      DO 70 I=1,N
      T=Z
      DO 60 J=1,M
   60 T=T+B(J,M)*A(I,J)
   70 X(I)=T
      CALL G(N,X,F)
      U=Z
      DO 80 I=1,N
   80 U=U+F(I)*F(I)
      IF (U.LT.EP1) GOTO 180
      S=Z
      DO 100 I=1,M
      T=B(I,M)
      DO 90 J=1,N
   90 T=T+B(I,J)*F(J)
      A(M,I)=T
      T=DABS(T)
      IF (T.GT.S) THEN
      S=T
      K=I
      ENDIF
  100 CONTINUE
      T=1/A(M,K)
      DO 130 J=1,M
      S=B(K,J)*T
      B(K,J)=S
      DO 110 I=K-1,1,-1
  110 B(I,J)=B(I,J)-A(M,I)*S
      DO 120 I=K+1,M
  120 B(I,J)=B(I,J)-A(M,I)*S
  130 CONTINUE
      IF (U.LT.P) THEN
      P=U
      L=K
      ENDIF
      T=Z
      DO 140 I=1,N
      S=A(I,K1)-X(I)
  140 T=T+S*S
      DO 150 I=1,N
  150 A(I,K)=X(I)
      K1=K
      IF (T.GT.EP) GOTO 50
      IF=1
  160 DO 170 I=1,N
  170 X(I)=A(I,L)
  180 CONTINUE
      CALL G(N,X,F)
      RETURN
      END
