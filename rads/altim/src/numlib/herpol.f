C    ****************************************************************** 
C                                                      15. HERPOL
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE VALUES OF HERMITE POLYNOMIALS HI(X)
C                                                                       
C    USAGE                                                              
C    CALL HERPOL(N1,X,H)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N1   I OPTION FOR POLYNOMIALS TO BE CALCULATED:
C           H(I-1) FOR I=1(1)N1, WHERE N1=N+1
C    X    I ARGUMENT VALUE
C    H    O ARRAY H(N1) CONTAINS THE RESULTANT VALUES HI
C           WITH H0 IN H(1) ETC.
C
C    REMARKS
C    1. FOR N<0, L0(X) IS CALCULATED
C    2. HERMITE POLYNOMIALS ARE USUALLY COMPUTED ON THE
C       INTERVAL (-1,+1) ALTHOUGH THE PROCEDURE IS NOT
C       RESTRICTED TO THIS RANGE, THERE IS DECREASING ACCURACY
C       FOR GREATER VALUES OF ABS(X).
C
C    *******************************************************************
      SUBROUTINE HERPOL(N1,X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(N1)
      DATA U,T/1.D0,2.D0/
      PI=DSQRT(3.141592653589793D0)
      H(1)=U/DSQRT(PI)
      IF (N1.LE.1) GOTO 40
      H(2)=X*H(1)*DSQRT(T)
      S=U+X
      DO 20 I=3,N1
  20  H(I)=X*H(I-1)/DSQRT((I-1)/T)-DSQRT(U-U/(I-1))*H(I-2)
  40  CONTINUE
      RETURN
      END
