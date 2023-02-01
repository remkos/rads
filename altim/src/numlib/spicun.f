C    ****************************************************************** 
C                                                      65. SPICUN
C                                                                       
C    PURPOSE                                                            
C    TO INTERPOLATE, DIFFERENTIATE AND INTEGRATE A TABULATED
C    FUNCTION ON NON-EQUIDISTAND POINTS BY CUBIC SPLINES
C                                                                       
C    USAGE                                                              
C    CALL SPICUN(N,M,G,A,X,Z,F,F1,V1,VN,C,D)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I M=-1 GIVEN POINTS ARE X1....XN, N>5
C         I M>=1 UPPER BOUNDERY OF QUADRATURE IS XN
C    M    I M=-1 ARRAY A IS FILLED, MUST OCCUR ONCE
C         I M<-1 F(X) AND F'(X) ARE CALCULATED
C         I M>=1 INTEGRATION ON INTERVAL (XM,XN)
C    G    I ARRAY G(N) CONTAINS THE FUNCTIONVALUES
C    A    O ARRAY A(N) CONTAINS THE MOMENTS AT POINTS XK
C           AND IS FILLED IF (M.EQ.-1), AND IS USED WHEN (M.NE.-1)
C    X    I X(N) CONTAINS THE ABCISSAE FOR WHICH F(X) AND F'(X)
C           ARE CALCULATED IF M<-1
C    Z    I IF M<-1 F(Z) AND F'(Z) ARE CALCULATED
C    F    O M<-1 VALUE F(Z)
C         O M>=1 VALUE OF QUADRATURE
C    F1   O M<-1 VALUE F'(Z)
C    V1     V1(N) AUXILIARY WORKINGSPACE
C    VN     VN(N) AUXILIARY WORKINGSPACE
C    C      C(N)  AUXILIARY WORKINGSPACE
C    D      D(N)  AUXILIARY WORKINGSPACE
C
C    SUBROUTINES REQUIRED
C    SPICUN USES SUBROUTINE SPICU1
C
C    REMARKS
C    IF SPICUN IS USED, THEN FIRST THE ARRAY A MUST BE FILLED
C    WITH THE MOMENTS. THIS IS DONE BY EXECUTING SPICUN WITH
C    M=-1 ONCE. AFTER THIS SPICUN CAN BE USED FOR INTERPOLATION
C    DIFFERENTIATION AND INTERGRATION.
C
C    METHOD
C    REFER TO RC-TWA-76001 OF COMPUTER CENTRE D.U.T.
C
C    *******************************************************************
      SUBROUTINE SPICUN(N,M,G,A,X,Z,F,F1,V1,VN,C,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(N),A(N),X(N),V1(N),VN(N),C(N),D(N)
      DATA ZZ,UU/0.D0,1.D0/
      N1=N-1
      IF (M.EQ.-1) THEN
      DO 10 I=2,N
   10 A(I)=6.D0*(G(I)-G(I-1))/(X(I)-X(I-1))
      DO 20 I=2,N1
   20 A(I)=A(I+1)-A(I)
      T=2.D0*(X(3)-X(1))
      D(2)=UU/T
      DO 30 I=2,N-2
      C(I)=(X(I+1)-X(I))*D(I)
      T=2.D0*(X(I+2)-X(I))-T*C(I)*C(I)
   30 D(I+1)=UU/T
      DO 40 I=2,N1
      V1(I)=ZZ
   40 VN(I)=ZZ
      V1(2)=UU
      VN(N1)=UU
      CALL SPICU1(N,1,V1,C,D)
      CALL SPICU1(N,0,VN,C,D)
      CALL SPICU1(N,1,A,C,D)
      R=ZZ
      T=ZZ
      S=ZZ
      V=ZZ
      P=ZZ
      Q=ZZ
      DO 70 I=2,5
      U=UU
      W=UU
      DO 50 J=2,5
      IF (J.NE.I) THEN
      U=U*(X(1)-X(J))
      W=W*(X(I)-X(J))
      ENDIF
   50 CONTINUE
      U=U/W
      R=R+U*V1(I)
      V=V+U*VN(I)
      P=P+U*A(I)
      U=UU
      W=UU
      DO 60 J=1,4
      IF (J.NE.I-1) THEN
      U=U*(X(N)-X(N-J))
      W=W*(X(N-I+1)-X(N-J))
      ENDIF
   60 CONTINUE
      U=U/W
      S=S+U*V1(N-I+1)
      T=T+U*VN(N-I+1)
      Q=Q+U*A(N-I+1)
   70 CONTINUE
      W=X(2)-X(1)
      U=X(N)-X(N1)
      R=UU+W*R
      V=U*V
      S=W*S
      T=UU+U*T
      E=R*T-S*V
      A(1)=(T*P-Q*V)/E
      A(N)=(R*Q-P*S)/E
      R=W*A(1)
      S=U*A(N)
      DO 80 I=2,N1
   80 A(I)=A(I)-R*V1(I)-S*VN(I)
      ELSE
      IF (M.LT.-1) THEN
      DO 90 I=N,1,-1
      K=I
      IF (Z.GE.X(K)) GOTO 100
   90 CONTINUE
  100 I=K+1
      W=X(I)-X(K)
      T=(Z-X(K))/W
      S=T*(UU-T)*W*W/6.D0
      R=A(I)*S
      S=A(K)*S
      F=T*(G(I)-R)+(UU-T)*(G(K)-S)-(R+S)
      S=T*T*A(I)-(UU-T)*(UU-T)*A(K)-(A(I)-A(K))/3.D0
      F1=(G(I)-G(K))/W+0.5D0*W*S
      ELSE
      S=ZZ
      T=ZZ
      DO 110 I=M+1,N
      K=I-1
      W=X(I)-X(K)
      S=S+(A(I)+A(K))*W*W*W
  110 T=T+(G(I)+G(K))*W
      F=0.5D0*T-S/24.D0
      ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE SPICU1(N,L,B,C,D)
      DOUBLE PRECISION B,C,D
      DIMENSION B(N),C(N),D(N)
      N1=N-1
      IF (L.EQ.1) THEN
      DO 10 I=3,N1
   10 B(I)=B(I)-C(I-1)*B(I-1)
      DO 20 I=2,N1
   20 B(I)=B(I)*D(I)
      ELSE
      B(N1)=D(N1)
      ENDIF
      DO 30 I=N-2,2,-1
   30 B(I)=B(I)-C(I)*B(I+1)
      RETURN
      END
