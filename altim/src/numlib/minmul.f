C    ****************************************************************** 
C                                                      52. MINMUL
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE ROOT OF A REAL FUNCTION FUN(X) ON AN
C    INTERVAL (A,B), WHERE FUN(A)*FUN(B)<0
C                                                                       
C    USAGE                                                              
C    VAR=MINMUL(A,B,EP,FUN)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I THE LOWER BOUNDARY OF THE INTERVAL
C    B    I THE UPPER BOUNDARY OF THE INTERVAL
C    EP   I TOLERANCE VALUE, THE PROCESS IS ENDED
C           WHEN AN INTERVAL (X0,X1) IS OBTAINED WITH
C           FUN(X0)*FUN(X1)<=0 AND X1-X0<EP
C    FUN  I NAME OF THE GIVEN FUNCTION
C
C    REMARKS
C    THE METHOD USED IS BASED ON A RATIONAL APPROXIMATION WITH
C    SOME ADDITIONS TO INCREASE THE SPEED OF CONVERGENCE
C
C    REFERENCES
C    RC-TWA-75003 NUMLIBDA METHODEBESCHRIJVINGEN
C
C    *******************************************************************
      FUNCTION MINMUL(A,B,EP,FUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MINMUL
      EXTERNAL FUN
      DIMENSION F(2),X(2)
      DATA Z,U,D/0.D0,1.D0,0.5D0/
      F(1)=FUN(A)
      IF (F(1).EQ.Z) THEN
      S=A
      GOTO 100
      ENDIF
      F(2)=FUN(B)
      IF (F(2).EQ.Z) THEN
      S=B
      GOTO 100
      ENDIF
      X(1)=A
      X(2)=B
      R=B-A
    4 P=D
      Q=Z
   10 S=X(1)+P*R
      FS=FUN(S)
      IF (FS.EQ.Z) GOTO 100
      T=FS*F(2)
      IF (T.EQ.Z) THEN
      X(1)=S
      GOTO 90
      ENDIF
      J=0
      IF (T.GT.Z) J=1
      K=J+1
      T=X(K)-S
      X(K)=S
      S=F(K)
      F(K)=FS
      FS=S
      R=X(2)-X(1)
      IF (R.LT.EP) GOTO 90
      Q=Q+J-D
      IF (DABS(Q).GT.1.9) GOTO 4
      IF (DABS(T).GT.EP) THEN
      S=F(1)*(F(2)-FS)
      P=S/(S+F(2)*(FS-F(1))*(1-P*J)/(J+P*(1-J)))
      IF ((P.GT.U).OR.(P.LT.Z)) P=D
      ELSE
      P=0.01+0.98*J
      ENDIF
      GOTO 10
   90 S=D*(X(1)+X(2))
  100 MINMUL=S
      RETURN
      END
