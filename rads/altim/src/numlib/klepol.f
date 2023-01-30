C    ****************************************************************** 
C                                                      28. KLEPOL
C                                                                       
C    PURPOSE                                                            
C    TO APPROXIMATE A SET OF (X,Y) POINTS BY A POLYNOMIAL
C    OF WANTED DEGREE USING THE LEAST SQUARES METHOD
C                                                                       
C    USAGE                                                              
C    CALL KLEPOL(X,Y,M,K1,SIG,C)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    X    I X(M) CONTAINS THE ABSCISSA OF THE GIVEN POINTS
C    Y    I Y(M) CONTAINS THE ORDINATES OF THE GIVEN POINTS
C    M    I THE NUMBER OF POINTS
C    K1   I THE WANTED DEGREE OF THE POLYNOMIAL + 1
C    SIG  O SIG(K1), SEE UNDER REMARKS
C    P    O P(K1) CONTAINS THE COMPUTED COEFFICIENTS
C           OF THE POLYNOMIAL SUM (J:1,K1) P(J)*X**(J-1)
C    C      WORKING AREA OF DIMENSION OF AT LEAST 6+6*K+2*M
C
C    REMARKS
C    ARRAY SIG(K1) CONTAINS THE FOLLOWING QUANTITIES:
C     SIG(1)=D1/(M-1), SIG(I)=DI/(M-I-1) I=2(1)K1,
C    WHERE DI= SUM(J:1,M) (Y(J)-FI(X(J)))**2
C    AND WHERE FI(X) IS A POLYNOMIAL OF DEGREE I-1
C
C    REFERENCES
C    RC-TWA-75003 NUMLIBDA METHODEBESCHRIJVINGEN
C
C    *******************************************************************
      SUBROUTINE KLEPOL(X,Y,M,K1,SIG,P,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SWX
      DIMENSION X(M),Y(M),SIG(K1),P(K1),C(*)
      DATA Z,U/0.D0,1.D0/
      K=K1-1
      IC=K1+2
      IZ=IC+K1
      IB=IZ+K1
      IS=IB+K1
      IA=IS+K
      I1=IA+K
      I2=I1+M
      SWX=.TRUE.
      C(IB)=Z
      C(1)=Z
      C(2)=Z
      DEL=Z
      OME=Z
      C(IC)=U
      R2=M
      DO 10 I=1,M
      DEL=DEL+Y(I)*Y(I)
      C(I+I1)=Z
      C(I+I2)=U
   10 OME=OME+Y(I)
      C(IS)=OME/R2
      C(IZ)=C(IS)
      DEL=DEL-C(IS)*OME
      SIG(1)=DEL/(M-U)
      DO 80 I=0,K-1
      D=Z
      DO 20 J=1,M
   20 D=D+X(J)*C(J+I2)**2
      C(I+1+IA)=D/R2
      R1=R2
      R2=Z
      OME=Z
      RU=Z
      DO 30 J=1,M
      D=C(I+IB)*C(J+I1)
      C(J+I1)=C(J+I2)
      C(J+I2)=(X(J)-C(I+1+IA))*C(J+I2)-D
      R2=R2+C(J+I2)**2
      OME=OME+Y(J)*C(J+I2)
   30 RU=RU+X(J)*C(J+I2)*C(J+I1)
      C(I+1+IB)=RU/R1
      C(I+1+IS)=OME/R2
      DEL=DEL-C(I+1+IS)*OME
      SIG(I+2)=DEL/(M-I-1)
      IF ((SIG(I+2).LT.SIG(I+1)).AND.(SWX)) GOTO 60
      IF (SWX) GOTO 40
      IF (SIG(I+2).LE.G) SWX=.TRUE.
      GOTO 60
   40 CONTINUE
      DO 50 J=0,I
   50 P(J+1)=C(J+IZ)
      SWX=.FALSE.
      G=SIG(I+1)
      N=I
   60 CONTINUE
      DO 70 J=0,I
      D=C(J+2)*C(I+IB)
      C(J+2)=C(J+IC)
      C(J+IC)=C(J+1)-C(I+1+IA)*C(J+IC)-D
   70 C(J+IZ)=C(J+IZ)+C(I+1+IS)*C(J+IC)
      C(I+1+IZ)=C(I+1+IS)
      C(I+1+IC)=U
      C(I+3)=Z
   80 CONTINUE
      IF (SWX) THEN
      DO 90 J=0,K
   90 P(J+1)=C(J+IZ)
      N=K
      ENDIF
      DO 100 I=N+1,K
  100 P(I+1)=Z
      RETURN
      END
