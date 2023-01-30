C    ****************************************************************** 
C                                                      11. EIGJAD
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE N EIGENVALUES OF A REAL MATRIX A,
C    WHERE A IS A SYMMETRIC MATRIX IN VECTOR FORM
C                                                                       
C    USAGE                                                              
C    CALL EIGJAD(A,N,RHO,H)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX A. THE DIMENSION IS AT LEAST N*(N+1)/2
C           A CAN BE FORMED BY SUBROUTINE MATSYD
C         O THE EIGENVALUES ARE GENERATED ON THE MAIN
C           DIAGONAL A(I+H(I))
C    N    I DIMENSION OF THE ORIGINAL MATRIX
C    RHO  I TOLERANCE VALUE
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH MATSY1
C                                                                       
C    METHOD                                                             
C    JACOBI ITERATION
C                                                                       
C    *******************************************************************
      SUBROUTINE EIGJAD(A,N,RHO,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION NORM1,NORM2,MU,INT1
      INTEGER P,Q,H
      DIMENSION A(*),H(N)
      DATA Z,U/0.0D0,1.0D0/
      INT1=Z
      DO 30 I=2,N
      DO 30 J=1,I-1
      K=J+H(I)
   30 INT1=INT1+A(K)*A(K)
      INT1=IN1+INT1
      NORM1=DSQRT(INT1)
      NORM2=(RHO/N)*NORM1
      THR=NORM1
      IND=0
   40 THR=THR/N
   50 CONTINUE
      DO 80 Q=2,N
      J=H(Q)
      DO 70 P=1,Q-1
      V2=A(P+J)
      IF (DABS(V2).LT.THR) GOTO 70
      IND=1
      V1=A(P+H(P))
      V3=A(Q+J)
      MU=0.5D0*(V1-V3)
      IF (MU.EQ.Z) THEN
      OMEGA=U
      ELSE
      OMEGA=DSIGN(U,MU)
      ENDIF
      OMEGA=-OMEGA*V2/DSQRT(V2*V2+MU*MU)
      SINT=U+DSQRT(U-OMEGA*OMEGA)
      SINT=OMEGA/DSQRT(SINT+SINT)
      MU=SINT*SINT
      OMEGA=1-MU
      COST=DSQRT(OMEGA)
      DO 54 I=1,N
      IF (I.LT.P) THEN
      K=I+H(P)
      ELSE
      K=P+H(I)
      ENDIF
      IF (I.LT.Q) THEN
      L=I+H(Q)
      ELSE
      L=Q+H(I)
      ENDIF
      INT1=A(K)*COST-A(L)*SINT
      A(L)=A(K)*SINT+A(L)*COST
   54 A(K)=INT1
      OMEGA=COST*COST
      INT1=COST*SINT
      A(P+H(P))=V1*OMEGA+V3*MU-2*V2*INT1
      A(Q+J)=V1*MU+V3*OMEGA+2*V2*INT1
      A(P+J)=(V1-V3)*INT1+V2*(OMEGA-MU)
   70 CONTINUE
   80 CONTINUE
      IF (IND.EQ.1) THEN
      IND=0
      GOTO 50
      ELSEIF (THR.GT.NORM2) THEN
      GOTO 40
      ENDIF
      RETURN
      END
