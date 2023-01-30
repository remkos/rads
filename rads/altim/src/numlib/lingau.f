C    ****************************************************************** 
C                                                      34. LINGAU
C                                                                       
C    PURPOSE                                                            
C    TO SOLVE A GENERAL SYSTEM OF LINEAR EQUATIONS AX=B                 
C                                                                       
C    USAGE                                                              
C    CALL LINGAU(N,L,A,B,D,EP,IF)                                       
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I NUMBER OF EQUATIONS                                         
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING             
C           (SUB)PROGRAM (L.GE.N)                                       
C    A    I REAL ARRAY OF DIMENSION (L,K) WHERE K.GE.N (DESTROYED)      
C    B    I RIGHT HAND SIDE, VECTOR OF DIMENSION AT LEAST N
C         O SOLUTION OF THE EQUATIONS                                   
C    D    O DETERMINANT OF A                                            
C    EP   I REAL TOLERANCE VALUE                                        
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT                         
C                              IF=1 PIVOT LESS EP
C                                                                       
C    REMARK                                                             
C    REALS ARE IN DOUBLE PRECISION
C                                                                       
C    METHOD                                                             
C    GAUSS ELIMINATION WITH COLUMN PIVOTTING                            
C                                                                       
C    *******************************************************************
      SUBROUTINE LINGAU(N,L,A,B,D,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(L,N),B(N)
      IF=0
      D=1.0D0
      DO 30 M=1,N
      P=0.0D0
      DO 8 I=M,N
      T=DABS(A(I,M))
      IF (P.GE.T) GOTO 8
      K=I
      P=T
    8 CONTINUE
      IF (P.GE.EP) GOTO 10
      IF=1
      GOTO 100
   10 IF (K.EQ.M) GOTO 14
      DO 12 J=M,N
      T=A(M,J)
      A(M,J)=A(K,J)
   12 A(K,J)=T
      T=B(M)
      B(M)=B(K)
      B(K)=T
      D=-D
   14 P=1/A(M,M)
      D=D*A(M,M)
      A(M,M)=P
      DO 20 I=M+1,N
      T=P*A(I,M)
      DO 16 J=M+1,N
   16 A(I,J)=A(I,J)-T*A(M,J)
   20 B(I)=B(I)-T*B(M)
   30 CONTINUE
      DO 36 I=N,1,-1
      T=B(I)
      DO 34 J=I+1,N
   34 T=T-A(I,J)*B(J)
   36 B(I)=T*A(I,I)
  100 RETURN
      END
