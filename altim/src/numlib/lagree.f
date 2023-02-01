C    ****************************************************************** 
C                                                      30. LAGREE
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE SUM OF LAGUERRE SERIES
C    A0*L0(X).........AN*LN(X), FOR GIVEN COEFFICIENTS AI
C                                                                       
C    USAGE                                                              
C    CALL LAGREE(N1,X,A,Y)
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
C
C    *******************************************************************
      SUBROUTINE LAGREE(N1,X,A,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N1)
      DATA U,T/1.D0,2.D0/
      Y=0.0D0
      IF (N1.LT.1) GOTO 40
      S=U+X
      PB=Y
      PA=Y
      DO 20 I=N1,2,-1
      PC=(T-S/I)*PB+A(I)-(U-U/(I+1))*PA
      PA=PB
  20  PB=PC
      Y=PB*(U-X)-PA/T+A(1)
  40  CONTINUE
      RETURN
      END
