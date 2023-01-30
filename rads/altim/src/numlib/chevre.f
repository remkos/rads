C    ****************************************************************** 
C                                                       7. CHEVRE
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE SUM OF SHIFTED CHEBYSHEV SERIES
C    A0*T$0(X).........AN*T$N(X), FOR GIVEN COEFFICIENTS AI
C                                                                       
C    USAGE                                                              
C    CALL CHEVRE(N1,X,A,Y)
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
C    2. SHIFTED CHEBYSHEV POLYNOMIALS ARE USUALLY COMPUTED ON THE
C       INTERVAL (0,+1) ALTHOUGH THE PROCEDURE IS NOT
C       RESTRICTED TO THIS RANGE, THERE IS DECREASING ACCURACY
C       FOR VALUES OF X OUTSIDE THIS INTERVAL.
C
C    *******************************************************************
      SUBROUTINE CHEVRE(N1,X,A,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N1)
      Y=0.0D0
      IF (N1.LT.1) GOTO 40
      U=X+X-1.D0
      PB=Y
      PA=Y
      UU=U+U
      DO 20 I=N1,2,-1
      PC=UU*PB+A(I)-PA
      PA=PB
  20  PB=PC
      Y=PB*U-PA+A(1)
  40  CONTINUE
      RETURN
      END
