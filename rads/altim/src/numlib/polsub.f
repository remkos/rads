C    ****************************************************************** 
C                                                      59. POLSUB
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE COEFFICIENTS OF A POLYNOMIAL C, WHICH IS
C    THE RESULT WHEN A POLYNOMIAL B IS USED AS THE ARGUMENT
C    OF ANOTHER POLYNOMIAL A.
C                                                                       
C    USAGE                                                              
C    CALL POLSUB(A,B,C,D,N,M)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I ARRAY A(N) CONTAINS THE COEFFICIENTS OF THE
C           POLYNOMIAL A: SUM A(J)*X**(J-1) J=1(1)N
C    B    I ARRAY B(M) CONTAINS THE COEFFICIENTS OF THE
C           POLYNOMIAL B: SUM B(J)*X**(J-1) J=1(1)M
C    C    O ARRAY C(K) CONTAINS THE COEFFICIENTS OF THE
C           RESULTANT POLYNOMIAL C: SUM C(J)*X**(J-1)
C           J=1(1)K, WHERE K=(N-1)*(M-1)+1
C    B    I VECTOR D IS A WORKING SPACE OF DIMENSION AT LEAST M*N
C    N    I DEGREE OF POLYNOMIAL A IS (N-1)
C    M    I DEGREE OF POLYNOMIAL B IS (M-1)
C
C    REMARKS
C    POLSUB USES SUBROUTINE POLVER
C
C    *******************************************************************
      SUBROUTINE POLSUB(A,B,C,D,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(M),C(*),D(*)
      D(1)=A(N)
      K=1
      DO 20 I=N-1,1,-1
      CALL POLVER(D,B,C,K,M)
      C(1)=C(1)+A(I)
      K=K+M-1
      DO 10 J=1,K
   10 D(J)=C(J)
   20 CONTINUE
      RETURN
      END
