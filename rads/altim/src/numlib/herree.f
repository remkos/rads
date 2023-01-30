C    ****************************************************************** 
C                                                      16. HERREE
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE SUM OF HERMITE SERIES
C    A0*H0(X).........AN*HN(X), FOR GIVEN COEFFICIENTS AI
C                                                                       
C    USAGE                                                              
C    CALL HERREE(N1,X,A,Y)
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
      SUBROUTINE HERREE(N1,X,A,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N1)
      DATA U,T,PI/1.D0,2.D0,3.141592653589793D0/
      Y=0.0D0
      IF (N1.LT.1) GOTO 40
      P=DSQRT(PI)
      PB=Y
      PA=Y
      DO 20 I=N1,2,-1
      PC=X*PB/DSQRT(I/T)+A(I)-PA*DSQRT(I/(I+U))
      PA=PB
  20  PB=PC
      Y=((T*PB*X-PA)/DSQRT(T)+A(1))/DSQRT(P)
  40  CONTINUE
      RETURN
      END
