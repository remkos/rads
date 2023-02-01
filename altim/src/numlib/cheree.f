C    ****************************************************************** 
C                                                       3. CHEREE
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE SUM OF CHEBYSHEV SERIES
C    A0*T0(X).........AN*TN(X), FOR GIVEN COEFFICIENTS AI
C                                                                       
C    USAGE                                                              
C    CALL CHEREE(N1,X,A,Y)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N1   I MAXIMUM INDEX +1 IN THE SUM AND ORDER +1
C           OF THE POLYNOMIAL. (N1=N+1)
C    X    I ARGUMENT VALUE
C    A    I ARRAY A(N1) CONTAINS THE COEFFICIENTS AI
C    Y    O VALUE OF THE SUM
C
C    REMARKS
C    1. FOR N<0, Y=0 IS RETURNED
C    2. CHEBYSHEV POLYNOMIALS ARE USUALLY COMPUTED ON THE
C       INTERVAL (-1,+1) ALTHOUGH THE PROCEDURE IS NOT
C       RESTRICTED TO THIS RANGE, THERE IS DECREASING ACCURACY
C       FOR VALUES OF X OUTSIDE THIS INTERVAL.
C
C    *******************************************************************
      SUBROUTINE CHEREE(N1,X,A,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N1)
      Y=0.0D0
      IF (N1.LT.1) GOTO 40
      U=X+X
      PB=Y
      PA=Y
      DO 20 I=N1,2,-1
      PC=U*PB+A(I)-PA
      PA=PB
  20  PB=PC
      Y=PB*X-PA+A(1)
  40  CONTINUE
      RETURN
      END
