C    ****************************************************************** 
C                                                      60. POLVER
C                                                                       
C    PURPOSE                                                            
C    TO MULTIPLICATE TWO POLYNOMIALS A*B=C
C                                                                       
C    USAGE                                                              
C    CALL POLVER(A,B,C,N,M)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I ARRAY A(N) CONTAINS THE COEFFICIENTS OF THE
C           POLYNOMIAL A: SUM A(J)*X**(J-1) J=1(1)N
C    B    I ARRAY B(M) CONTAINS THE COEFFICIENTS OF THE
C           POLYNOMIAL B: SUM B(J)*X**(J-1) J=1(1)M
C    C    O ARRAY C(K) CONTAINS THE COEFFICIENTS OF THE
C           POLYNOMIAL C: SUM C(J)*X**(J-1) J=1(1)N+M
C    N    I DEGREE OF POLYNOMIAL A IS (N-1)
C    M    I DEGREE OF POLYNOMIAL B IS (M-1)
C
C    REMARKS
C    DEGREE K OF POLYNOMIAL C IS (N+M-2)
C
C    *******************************************************************
      SUBROUTINE POLVER(A,B,C,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(M),C(*)
      DO 10 I=1,M+N-1
   10 C(I)=0.D0
      DO 30 J=1,M
      DO 20 I=J,N+J-1
   20 C(I)=C(I)+B(J)*A(I-J+1)
   30 CONTINUE
      RETURN
      END
