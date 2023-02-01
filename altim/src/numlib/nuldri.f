C    ****************************************************************** 
C                                                      54. NULDRI
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE ROOTS OF A CUBIC EQUATION
C                                                                       
C    USAGE                                                              
C    CALL NULDRI(A,X,IE)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I ARRAY A(3) CONTAINS THE COEFFICIENTS OF THE CUBIC
C           EQUATION X*X*X+A(3)*X*X+A(2)*X+A(1)
C    X    I ARRAY X(3) CONTAINS THE ROOTS, SEE IE
C    IE   O EXITCODE
C           IE=0: REAL ROOTS, WHERE X(1)>=X(2)>=X(3)
C           IE=1: REAL ROOT X(1) AND COMPLEX ROOTS
C                 X(2)+IX(3) AND X(2)-IX(3)
C           IE=2: THE RESULT IS QUESTIONABLE
C
C    REMARKS
C    NULDRI USES SUBROUTINE NULTWE
C
C    REFERENCES
C    RC-TWA-81002 ROOTS OF QUADRATIC AND CUBIC EQUATIONS
C
C    *******************************************************************
      SUBROUTINE NULDRI(A,X,IE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3),X(3)
      DATA Z,U,T,D,E/0.D0,1.D0,2.D0,3.D0,1.D-15/
      Y=Z
      IF (A(1).EQ.Z) GOTO 20
      P=A(3)*A(3)-D*A(2)
      Y=-A(3)/D
      IF (P.GT.Z) THEN
      R=A(2)/D
      CALL NULTWE(Y,R,IE)
      W=((R+A(3))*R+A(2))*R+A(1)
      S=((Y+A(3))*Y+A(2))*Y+A(1)
      IF (W.GT.-S) Y=R
      ENDIF
      S=((Y+A(3))*Y+A(2))*Y+A(1)
      IF (S.EQ.Z) GOTO 20
      P=DSIGN(U,S)*DABS(S)**(U/D)
      Y=Y-P
      W=E*DABS(A(1))
   10 P=(D*Y+T*A(3))*Y+A(2)
      IF (P.EQ.Z) THEN
      IE=2
      GOTO 20
      ENDIF
      S=((Y+A(3))*Y+A(2))*Y+A(1)
      IF (DABS(S).LE.W) GOTO 20
      P=S/P
      Y=Y-P
      IF (DABS(P).GT.E*(U+DABS(Y))) GOTO 10
   20 B=A(3)+Y
      IF (DABS(B).LT.0.2*DABS(Y)) THEN
      C=-A(1)/Y
      B=(A(2)-C)/(T*Y)
      ELSE
      B=-B/T
      C=A(2)-T*B*Y
      ENDIF
      CALL NULTWE(B,C,IE)
      X(1)=Y
      X(2)=B
      X(3)=C
      IF ((IE.EQ.0).AND.(Y.LT.B)) THEN
      X(1)=B
      X(2)=C
      X(3)=Y
      ENDIF
      RETURN
      END
