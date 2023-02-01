C    ****************************************************************** 
C                                                     49B. MATSYR
C                                                                       
C    PURPOSE                                                            
C    TO TRANSFORM A VECTORIZED MATRIX A INTO A SYMMETRIC MATRIX C
C    (THIS IS THE REVERSE OF MATSYD)
C                                                                       
C    USAGE                                                              
C    CALL MATSYR(N,L,A,C)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I DIMENSION OF A
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING             
C           (SUB)PROGRAM (L.GE.N)                                       
C    A    I VECTOR OF DIMENSION AT LEAST N*(N+1)/2
C    C    O RESULTANT REAL ARRAY OF DIMENSION (L,N)
C                                                                       
C    *******************************************************************
      SUBROUTINE MATSYR(N,L,A,C)
      REAL*8 A(*),C(L,N)
      K=0
      DO 20 J=1,N
         DO 20 I=1,J
            K=K+1
            C(J,I)=A(K)
   20       C(I,J)=A(K)
      END
