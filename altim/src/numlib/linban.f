C    ****************************************************************** 
C                                                      32. LINBAN
C                                                                       
C    PURPOSE                                                            
C    TO SOLVE A GENERAL SYSTEM OF LINEAR EQUATIONS AX=B,
C    WHERE A IS A BAND MATRIX
C                                                                       
C    USAGE                                                              
C    CALL LINBAN(M,N,K,L,A,B,D,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    M    I WIDTH OF THE BAND OF A
C    N    I NUMBER OF EQUATIONS                                         
C    K    I NUMBER OF LOWER CODIAGONALS
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING             
C           (SUB)PROGRAM (L.GE.N)                                       
C    A    I REAL ARRAY OF DIMENSION (L,M) SEE REMARKS. (DESTROYED)
C    B    I RIGHT HAND SIDE, REAL ARRAY OF DIMENSION AT LEAST N         
C         O SOLUTION OF THE EQUATIONS                                   
C    D    O DETERMINANT OF A                                            
C    EP   I REAL TOLERANCE VALUE                                        
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT                         
C                              IF=1 PIVOT LESS EP
C                                                                       
C    REMARKS
C    1. THE SIGNIFICANT ELEMENTS OF THE ORIGINAL MATRIX MUST BE
C       SHIFTED TO THE LEFT AS FAR AS POSSIBLE. THE EMPTY
C       POSITIONS MUST BE FILLED WITH ZEROES. THE DIAGONAL ELEMENTS
C       WILL THEN BE SITUATED IN THE FOLLOWING POSITIONS.
C          FOR  0<I<K+1 A(I,I)
C          FOR  K<I<N+1 A(I,K+1)
C    2. A CAN BE FILLED WITH THE AID OF SUBROUTINE MATBAN.
C                                                                       
C    METHOD                                                             
C    GAUSS ELIMINATION WITH COLUMN PIVOTTING                            
C                                                                       
C    *******************************************************************
      SUBROUTINE LINBAN(M,N,K,L,A,B,D,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER R,P,Q,U
      DIMENSION A(L,M),B(N)
      IF=0
      D=1.0D0
      R=M
      Q=N-M+1
      DO 30 I=1,N
      IF (I.GT.Q) R=R-1
      IF (K.NE.N) K=K+1
      T=0.0D0
      DO 8 U=I,K
      S=DABS(A(U,1))
      IF (T.GE.S) GOTO 8
      T=S
      P=U
    8 CONTINUE
      IF (T.GE.EP) GOTO 10
      IF=1
      GOTO 100
   10 IF (P.EQ.I) GOTO 14
      DO 12 J=1,R
      S=A(I,J)
      A(I,J)=A(P,J)
   12 A(P,J)=S
      S=B(I)
      B(I)=B(P)
      B(P)=S
      D=-D
   14 D=D*A(I,1)
      T=1/A(I,1)
      DO 16 J=2,R
   16 A(I,J)=T*A(I,J)
      B(I)=T*B(I)
      DO 24 U=I+1,K
      S=A(U,1)
      DO 20 J=2,R
   20 A(U,J-1)=A(U,J)-S*A(I,J)
      B(U)=B(U)-S*B(I)
   24 A(U,R)=0.0D0
   30 CONTINUE
      U=1
      DO 36 I=N-1,1,-1
      IF (U.LT.M) U=U+1
      DO 34 J=2,U
   34 B(I)=B(I)-A(I,J)*B(I-1+J)
   36 CONTINUE
  100 RETURN
      END
