C    ****************************************************************** 
C                                                       5. CHEVER
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE VALUES OF SHIFTED CHEBYSHEV POLYNOMIALS T$N(X)
C                                                                       
C    USAGE                                                              
C    CALL CHEVER(N1,X,T)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N1   I OPTION FOR POLYNOMIALS TO BE CALCULATED:
C           T$(I-1) FOR I=1(1)N1, WHERE N1=N+1
C    X    I ARGUMENT VALUE
C    T    O ARRAY T(N1) CONTAINS THE RESULTANT VALUES
C           WITH T0 IN T(1) ETC.
C
C    REMARKS
C    1. FOR N<0, THE SUBROUTINE HAS NO EFFECT
C    2. SHIFTED CHEBYSHEV POLYNOMIALS ARE USUALLY COMPUTED ON THE
C       INTERVAL (0,+1) ALTHOUGH THE PROCEDURE IS NOT
C       RESTRICTED TO THIS RANGE, THERE IS DECREASING ACCURACY
C       FOR VALUES OF X OUTSIDE THIS INTERVAL.
C
C    *******************************************************************
      SUBROUTINE CHEVER(N1,X,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N1)
      T(1)=1.0D0
      IF (N1.LE.1) GOTO 40
      T(2)=X*2.D0-1.D0
      U=T(2)+T(2)
      DO 20 I=3,N1
  20  T(I)=U*T(I-1)-T(I-2)
  40  CONTINUE
      RETURN
      END
