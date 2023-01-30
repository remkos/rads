C    ****************************************************************** 
C                                                      10. EIGJAC
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE N EIGENVALUES AND EIGENVECTORS OF
C    A REAL SYMMETRIC MATRIX A
C                                                                       
C    USAGE                                                              
C    CALL EIGJAC(A,S,N,L,RHO)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I THE N GIVEN REAL SYMMETRIC MATRIX, DEMENSION (N,N)
C         O THE EIGENVALUES ARE GENERATED ON THE MAIN DIAGONAL A(I,I)
C    S    O THE N RESULTANT EIGENVECTORS, DIMENSION (N,N)
C    N    I ORDER OF THE MATRIX A
C    L    I FIRST DIMENSION OF A AND S AS DECLARED IN THE CALLING
C           (SUB)PROGRAM. LA IS AT LEAST N
C    RHO  I TOLERANCE VALUE
C                                                                       
C    METHOD                                                             
C    JACOBI ITERATION
C                                                                       
C    *******************************************************************
      SUBROUTINE EIGJAC(A,S,N,L,RHO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION NORM1,NORM2,MU,INT1,IN
      INTEGER P,Q
      DIMENSION A(L,N),S(L,N)
      DATA Z,U/0.0D0,1.0D0/
      DO 20 I=1,N
      DO 10 J=1,N
   10 S(I,J)=Z
   20 S(I,I)=U
      INT1=Z
      DO 30 I=2,N
      DO 30 J=1,I-1
   30 INT1=INT1+2.D0*A(I,J)*A(I,J)
      NORM1=DSQRT(INT1)
      NORM2=(RHO/N)*NORM1
      THR=NORM1
      IND=0
   40 THR=THR/N
   50 CONTINUE
      DO 80 Q=2,N
      DO 70 P=1,Q-1
      IF (DABS(A(P,Q)).LT.THR) GOTO 70
      IND=1
      V1=A(P,P)
      V2=A(P,Q)
      V3=A(Q,Q)
      MU=0.5D0*(V1-V3)
      IF (DABS(MU).LT.1.0D-7*DABS(V2)) THEN
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
      INT1=A(I,P)*COST-A(I,Q)*SINT
      A(I,Q)=A(I,P)*SINT+A(I,Q)*COST
      A(I,P)=INT1
      INT1=S(I,P)*COST-S(I,Q)*SINT
      S(I,Q)=S(I,P)*SINT+S(I,Q)*COST
   54 S(I,P)=INT1
      DO 60 I=1,N
      A(P,I)=A(I,P)
   60 A(Q,I)=A(I,Q)
      INT1=COST*SINT
      IN=V2*INT1
      IN=IN+IN
      A(P,P)=V1*OMEGA+V3*MU-IN
      A(Q,Q)=V1*MU+V3*OMEGA+IN
      A(Q,P)=(V1-V3)*INT1+V2*(OMEGA-MU)
      A(P,Q)=A(Q,P)
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
