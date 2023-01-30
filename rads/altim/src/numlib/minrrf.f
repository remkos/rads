C    ****************************************************************** 
C                                                      53. MINRRF
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE ROOT OF A REAL FUNCTION FUN(X) ON AN
C    INTERVAL (A,B), WHERE FUN(A)*FUN(B)<0
C                                                                       
C    USAGE                                                              
C    CALL MINRRF(A,B,X,EP,MS,IF,FUN)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I THE LOWER BOUNDARY OF THE INTERVAL
C         O IF IF=2, THE LOWER BOUNDARY OF THE REMAINING INTERVAL
C    B    I THE UPPER BOUNDARY OF THE INTERVAL
C         O IF IF=2, THE UPPER BOUNDARY OF THE REMAINING INTERVAL
C    X    O RESULTANT ROOT OF THE FUNCTION FUN(X) ON (A,B)
C    EP   I TOLERANCE VALUE, THE PROCESS IS ENDED
C           WHEN AN INTERVAL (X0,X1) IS OBTAINED WITH
C           FUN(X0)*FUN(X1)<=0 AND X1-X0<EP
C    MS   I MAXIMUM NUMBER OF ITERATION STEPS
C    IF   O INTEGER ERRORCODE
C           IF=0 NORMAL EXIT, X IS THE MIDPOINT OF AN INTERVAL
C                WITH LENGTH LESS THAN EP
C           IF=1 NORMAL EXIT, X IS A POINT ,WHERE FUN(X)=0
C           IF=2 NO CONVERGENCE AFTER MS ITERATION STEPS.
C                IT IS POSSIBLE TO PROCEED BY CALLING MINRRF AGAIN
C           IF=3 SIGN(FUN(A)*FUN(B))=1
C           IF=4 SIGN(B-A)=-1
C    FUN  I NAME OF THE EXTERNAL FUNCTION USED
C
C    REMARKS
C    THE METHOD USED IS THE REVERSED REGULA FALSI.
C    IT IS VERY FAST FOR SIMPLE ROOTS.
C    MINRRF USES THE SUBROUTINES MINBER AND MINSUB
C
C    REFERENCES
C    RC-TWA-80009 REVERSED REGULA FALSI, COMPUTER CENTRE DELFT
C
C    *******************************************************************
      SUBROUTINE MINRRF(A,B,X,EP,MS,IF,FUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(4),Y(4)
      DATA U,H/1.D0,0.5D0/
      IF (B.LT.A) THEN
      IF=4
      GOTO 200
      ENDIF
      Y(1)=A
      Y(3)=B
      CALL MINBER(3,IF,X,F,Y,FUN)
      IF (IF.EQ.1) GOTO 200
      CALL MINBER(1,IF,X,F,Y,FUN)
      IF (IF.EQ.1) GOTO 200
      IF (F(3)*F(1).GT.0.) THEN
      IF=3
      GOTO 200
      ENDIF
      IF=0
      DO 120 M=1,MS
      D=Y(3)-Y(1)
      IF (D.LT.EP) THEN
      X=(Y(3)+Y(1))*H
      GOTO 200
      ENDIF
      Y(2)=(Y(3)*F(3)-Y(1)*F(1))/(F(3)-F(1))
      CALL MINBER(2,IF,X,F,Y,FUN)
      IF (IF.EQ.1) GOTO 200
      J=DSIGN(U,F(2))*DSIGN(U,F(3))
      T=Y(2)-F(2)*D/(F(3)-F(1))
      IF ((J*T.LE.J*Y(2-J)).OR.(J*T.GE.J*Y(2))) T=(Y(2)+Y(2-J))*H
      Y(4)=T
      CALL MINBER(4,IF,X,F,Y,FUN)
      IF (IF.EQ.1) GOTO 200
      I=DSIGN(U,F(4))*DSIGN(U,F(3))
      IF (I.EQ.J) THEN
      CALL MINSUB(J+2,4,F,Y)
      ELSE
      CALL MINSUB(J+2,2,F,Y)
      CALL MINSUB(2-J,4,F,Y)
      ENDIF
      D=(Y(3)-Y(1))
      S=-D*F(2+I)/F(2-I)
      IF (S.GE.D) GOTO 120
      IF (S.LT.H*EP) S=H*EP
      T=S
      DO 100 J=1,3
      IF (T.LT.D) THEN
      Y(2)=Y(2+I)-S*I
      CALL MINBER(2,IF,X,F,Y,FUN)
      IF (IF.EQ.1) GOTO 200
      IF (DSIGN(U,F(2))*DSIGN(U,F(2+I)).LT.0.D0) THEN
      CALL MINSUB(2-I,2,F,Y)
      GOTO 120
      ENDIF
      CALL MINSUB(2+I,2,F,Y)
      S=H*D/(4-J)
      T=T+S
      ENDIF
  100 CONTINUE
  120 CONTINUE
      IF=2
      A=Y(1)
      B=Y(3)
  200 CONTINUE
      RETURN
      END
      SUBROUTINE MINBER(I,IF,X,F,Y,FUN)
      DOUBLE PRECISION X,F,Y,FUN
      DIMENSION F(4),Y(4)
      EXTERNAL FUN
      F(I)=FUN(Y(I))
      IF (F(I).EQ.0.D0) THEN
      X=Y(I)
      IF=1
      ENDIF
      RETURN
      END
      SUBROUTINE MINSUB(I,J,F,Y)
      DOUBLE PRECISION F,Y
      DIMENSION F(4),Y(4)
      F(I)=F(J)
      Y(I)=Y(J)
      RETURN
      END
