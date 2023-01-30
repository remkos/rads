C    ****************************************************************** 
C                                                      49. MATSYD
C                                                                       
C    PURPOSE                                                            
C    TO TRANSFORM A SYMMETRIC MATRIX A INTO A VECTOR C in so-called
C    Upper Packed Storage Mode. Only the upper right triangle of
C    A has to be filled.
C    SUITED TO BE USED IN SUBROUTINE LINSYD
C                                                                       
C    USAGE                                                              
C    CALL MATSYD(N,L,A,C)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I DIMENSION OF A
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING             
C           (SUB)PROGRAM (L.GE.N)                                       
C    A    I REAL ARRAY OF DIMENSION (L,N)
C    C    O RESULTANT VECTOR OF DIMENSION AT LEAST N*(N+1)/2
C                                                                       
C    *******************************************************************
      SUBROUTINE MATSYD(N,L,A,C)
      REAL*8 A(L,N),C(*)
      K=0
      DO 20 J=1,N
         DO 20 I=1,J
            K=K+1
   20       C(K)=A(I,J)
      END
