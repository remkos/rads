C    ****************************************************************** 
C                                                      42. MATDET
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE DETERMINANT OF A MATRIX A
C                                                                       
C    USAGE                                                              
C    CALL MATDET(N,L,A,D,EP,IF)                                       
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I NUMBER OF EQUATIONS                                         
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING             
C           (SUB)PROGRAM (L.GE.N)                                       
C    A    I REAL ARRAY OF DIMENSION (L,K) WHERE K.GE.N (DESTROYED)      
C    D    O DETERMINANT OF A                                            
C    EP   I REAL TOLERANCE VALUE                                        
C    IF   O INTEGER ERROR CODE IF=0 NORMAL EXIT                         
C                              IF=1 PIVOT LESS EP
C                                                                       
C    METHOD                                                             
C    GAUSS ELIMINATION WITH COLUMN PIVOTTING                            
C                                                                       
C    *******************************************************************
      SUBROUTINE MATDET(N,L,A,D,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(L,N)
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
      D=-D
   14 P=1/A(M,M)
      D=D*A(M,M)
      A(M,M)=P
      DO 20 I=M+1,N
      T=P*A(I,M)
      DO 16 J=M+1,N
   16 A(I,J)=A(I,J)-T*A(M,J)
   20 CONTINUE
   30 CONTINUE
  100 RETURN
      END
