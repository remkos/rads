C    ****************************************************************** 
C                                                      37. LINTRI
C                                                                       
C    PURPOSE                                                            
C    TO SOLVE A LINEAR SYSTEM OF EQUATIONS AX=B, WHERE
C    A IS A REAL TRIDIAGONAL MATRIX
C                                                                       
C    USAGE                                                              
C    CALL LINTRI(N,C,L,U,X,F,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I ORDER OF THE MATRIX A
C    C    I C(N) CONTAINS THE MAINDIAGONAL ELEMENTS
C    U    I U(N) CONTAINS THE UPPER DIAGONAL ELEMENTS
C    L    I L(N) CONTAINS THE LOWER DIAGONAL ELEMENTS
C    X    I X(N) CONTAINS THE RIGHT HAND SIDE
C         O X(N) CONTAINS THE SOLUTION OF THE SYSTEM
C    V      V(N) IS AN AUXILIARY INTEGER VECTOR
C    F    0 A MEASURE OF THE ACCURACY, SEE REMARKS
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT
C                              IF=1 PIVOT=0
C
C    REMARKS
C    1. F IS THE QUOTIENT OF MIN(DABS(PIVOT)) AND
C       MAX(DABS(CI)+DABS(LI)+DABS(UI)) AND IS OF THE ORDER
C       OF THE ACCURACY OF THE SOLUTION. IF F<10**-15, THEN
C       THE SYSTEM WILL PROBABLY BE SINGULAR.
C    2. NOTE THAT L(1)=U(N)=0
C    3. AS A PIVOT AT STEP I MAX(DABS(CI),DABS(LI),DABS(UI))
C       IS CHOSEN.
C
C    *******************************************************************
      SUBROUTINE LINTRI(N,C,L,U,X,V,F,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION L
      INTEGER V
      DIMENSION C(N),L(N),U(N),X(N),V(N)
      DATA Z/0.D0/
      IF=1
      M=N-1
      R=Z
      FF=Z
      F=1.D70
      DO 10 I=1,N
      T=DABS(C(I))+DABS(L(I))+DABS(U(I))
      IF (T.GT.R) R=T
   10 CONTINUE
      DO 20 I=1,M
      V(I)=I
      T=C(I)
      J=I+1
      IF (DABS(T).LT.0.1*DABS(U(I))) THEN
      V(I)=0
      T=C(I)/L(J)
      L(I)=U(J)
      C(I)=L(J)
      IF (I.LT.N-1) U(J)=-L(J)*T
      S=X(J)
      X(J)=X(I)-S*T
      X(I)=S
      S=C(J)
      C(J)=U(I)-S*T
      U(I)=S
      ELSE
      IF (T.EQ.Z) GOTO 100
      C(J)=C(J)-L(J)*U(I)/T
      X(J)=X(J)-L(J)*X(I)/T
      ENDIF
      IF (DABS(C(I)).LT.F) F=DABS(C(I))
   20 CONTINUE
      IF (DABS(C(N)).LT.F) F=DABS(C(N))
      IF (C(N).EQ.Z) GOTO 100
      X(N)=X(N)/C(N)
      DO 30 I=M,1,-1
      IF ((V(I).EQ.0).AND.(I.LT.N-1)) THEN
      T=L(I)*X(I+2)
      ELSE
      T=Z
      ENDIF
   30 X(I)=(X(I)-U(I)*X(I+1)-T)/C(I)
      IF=0
      FF=F
  100 F=FF/R
      RETURN
      END
