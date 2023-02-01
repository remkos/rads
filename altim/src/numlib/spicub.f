C    ****************************************************************** 
C                                                      64. SPICUB
C                                                                       
C    PURPOSE                                                            
C    TO INTERPOLATE, DIFFERENTIATE AND INTEGRATE A TABULATED
C    FUNCTION ON EQUIDISTAND POINTS BY CUBIC SPLINES
C                                                                       
C    USAGE                                                              
C    CALL SPICUB(N,M,G,A,X1,X,H,F,F1)
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
C    X1   I LOWER BOUNDARY VALUE OF INTERVAL
C    X    I IF M<-1 F(X) AND F'(X) ARE CALCULATED
C    F    O M=-1 APPROXIMATION OF MAXIMUM ERROR IN F(X)
C         O M<-1 VALUE F(Z)
C         O M>=1 VALUE OF QUADRATURE
C    F1   O M=-1 APPROXIMATION OF MAXIMUM ERROR IN F'(X)
C         O M<-1 VALUE F'(Z)
C         O M>-1 APPROXIMATION OF ERROR IN THE QUADRATURE
C
C    SUBROUTINES REQUIRED
C    NONE
C
C    REMARKS
C    IF SPICUB IS USED, THEN FIRST THE ARRAY A MUST BE FILLED
C    WITH THE MOMENTS. THIS IS DONE BY EXECUTING SPICUN WITH
C    M=-1 ONCE. AFTER THIS SPICUN CAN BE USED FOR INTERPOLATION
C    DIFFERENTIATION AND INTERGRATION.
C
C    METHOD
C    REFER TO RC-TWA-76001 OF COMPUTER CENTRE D.U.T.
C
C    *******************************************************************
      SUBROUTINE SPICUB(N,M,G,A,X1,X,H,F,F1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(N),A(N)
      DATA Z,E,B,D/0.D0,1.D0,0.5D0,2.D0/
      U=B/(DSQRT(D)+E)
      IF (M.EQ.-1) THEN
      T=Z
      R=(E+B*U)/6.D0
      DO 10 K=N-2,3,-1
      W=G(K-1)+G(K+1)
      S=G(K-2)+G(K+2)-4.D0*W+6.D0*G(K)
      IF (DABS(S).GT.T) T=DABS(S)
      S=S*R
      A(K)=W-D*G(K)-S
      IF (K.EQ.N-2) F=S
   10 CONTINUE
      A(2)=G(1)+G(3)-D*G(2)-S
      A(N-1)=G(N-2)+G(N)-D*G(N-1)-F
      A(1)=3.D0*G(1)-9.D0*G(2)+10.D0*G(3)-5.D0*G(4)+G(5)-S
      A(N)=3.D0*G(N)-9.D0*G(N-1)+10.D0*G(N-2)-5.D0*G(N-3)+G(N-4)-F
      F=U*U*T/96.D0
      F1=U*T/(24.D0*H)
      ELSE
      IF (M.LT.-1) THEN
      T=(X-X1)/H
      I=T
      T=T-I
      I=I+2
      IF (I.GT.N) THEN
      I=N
      T=E
      ENDIF
      K=I-1
      S=T*(E-T)/6.D0
      R=A(I)*S
      S=A(K)*S
      F=T*(G(I)-R)+(E-T)*(G(K)-S)-(R+S)
      F1=G(I)-G(K)-(A(I)-A(K))/6.D0
      F1=(F1+B*(T*T*A(I)-(E-T)*(E-T)*A(K)))/H
      ELSE
      S=Z
      T=Z
      K=N-1
      DO 30 I=M+1,K
      S=S+A(I)
   30 T=T+G(I)
      F=H*(T+B*(G(N)+G(M))-(S+B*(A(N)+A(M)))/12.D0)
      F1=2.5D0*(G(M)+G(N))-9.D0*(G(M+1)+G(N-1))
      F1=12.D0*(G(M+2)+G(N-2))-7.D0*(G(M+3)+G(N-3))+F1
      F1=(1.5D0*(G(M+4)+G(N-4))+F1)*H*(0.2D0-U)/144.D0
      ENDIF
      ENDIF
      RETURN
      END
