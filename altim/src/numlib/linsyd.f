C    ****************************************************************** 
C                                                      36. LINSYD
C                                                                       
C    PURPOSE                                                            
C    TO SOLVE A SYSTEM OF LINEAR EQUATIONS AX=B
C    WHERE A IS A SYMMETRIC MATRIX IN VECTOR FORM
C                                                                       
C    USAGE                                                              
C    CALL LINSYD(N,A,B,H,F,X,D,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I NUMBER OF EQUATIONS                                         
C    A    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX. THE DIMENSION IS AT LEAST N*(N+1)/2
C    B    I RIGHT HAND SIDE, VECTOR OF DIMENSION AT LEAST N
C         O SOLUTION OF THE EQUATIONS                                   
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH MATSY1
C    F      INTEGER VECTOR WORKING SPACE OF DIMENSION AT LEAST N
C    X      REAL VECTOR WORKING SPACE OF DIMENSION AT LEAST
C    D    O DETERMINANT OF A                                            
C    EP   I REAL TOLERANCE VALUE                                        
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT                         
C                              IF=1 PIVOT LESS EP
C                                                                       
C    SUBROUTINES USED
C    LINSYD USES SUBROUTINE LINSY1
C                                                                       
C    METHOD                                                             
C    GAUSS ELIMINATION WITH PIVOTTING ON THE MAIN DIAGONAL
C                                                                       
C    *******************************************************************
      SUBROUTINE LINSYD(N,A,B,H,F,X,D,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER H,S,F
      DIMENSION A(*),B(N),H(N),F(N),X(N)
      DATA Z,U/0.0D0,1.0D0/
      IF=0
      D=U
      DO 30 S=1,N
         Y=Z
         DO 4 J=S,N
            T=DABS(A(J+H(J)))
            IF (Y.GE.T) GOTO 4
            L=J
            Y=T
   4     CONTINUE
         IF (Y.LT.EP) THEN
            IF=1
            GOTO 100
         ENDIF
  10     F(S)=L
         IF (L.NE.S) THEN
            DO 12 J=1,N
  12           CALL LINSY1(S,L,J,N,A,H)
            CALL LINSY1(S,L,S,N,A,H)
            T=B(L)
            B(L)=B(S)
            B(S)=T
         ENDIF
  14     Y=A(S+H(S))
         D=D*Y
         Y=U/Y
         DO 15 J=S+1,N
            T=A(S+H(J))
            X(J)=T
  15        A(S+H(J))=T*Y
         B(S)=B(S)*Y
         DO 20 I=S+1,N
            T=X(I)
            DO 16 J=I,N
  16           A(I+H(J))=A(I+H(J))-T*A(S+H(J))
  20        B(I)=B(I)-T*B(S)
  30  CONTINUE
      DO 36 S=N,1,-1
         T=B(S)
         DO 34 J=S+1,N
  34        T=T-A(S+H(J))*B(J)
  36     B(S)=T
      DO 40 J=N,1,-1
         M=F(J)
         T=B(M)
         B(M)=B(J)
  40     B(J)=T
 100  END
