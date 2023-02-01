C    ****************************************************************** 
C                                                      58. POLHOR
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE VALUE OF A POLYNOMIAL WITH REAL COEFFICIENTS
C                                                                       
C    USAGE                                                              
C    VAR=POLHOR(N1,A,X)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N1   I DEGREE +1 OF THE POLYNOMIAL
C    A    I ARRAY A(N1) CONTAINS THE COEFFICIENTS
C           A(I-1) FOR I=1(1)N1, WHERE N1=N+1
C    X    I ARGUMENT VALUE
C
C    REMARKS
C    1. POLHOR IS A DOUBLE PRECISION REAL FUNCTION
C    2. THE METHOD USEDIS HORNERS SCHEME.
C
C    *******************************************************************
      FUNCTION POLHOR (N1,A,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N1)
      T=A(N1)
      N=N1-1
      DO 20 I=N,1,-1
  20  T=X*T+A(I)
      POLHOR=T
      RETURN
      END
