C    ****************************************************************** 
C                                                      29. LAGPOL
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE VALUES OF LAGUERRE POLYNOMIALS LI(X)
C                                                                       
C    USAGE                                                              
C    CALL LAGPOL(N1,X,H)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N1   I OPTION FOR POLYNOMIALS TO BE CALCULATED:
C           H(I-1) FOR I=1(1)N1, WHERE N1=N+1
C    X    I ARGUMENT VALUE
C    H    O ARRAY H(N1) CONTAINS THE RESULTANT VALUES LI
C           WITH L0 IN H(1) ETC.
C
C    REMARKS
C    1. FOR N<0, L0(X) IS CALCULATED
C    2. LAGUERRE POLYNOMIALS ARE USUALLY COMPUTED ON THE
C       INTERVAL (-1,+1) ALTHOUGH THE PROCEDURE IS NOT
C       RESTRICTED TO THIS RANGE, THERE IS DECREASING ACCURACY
C       FOR GREATER VALUES OF ABS(X).
C
C    *******************************************************************
      SUBROUTINE LAGPOL(N1,X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(N1)
      DATA U,T/1.D0,2.D0/
      H(1)=U
      IF (N1.LE.1) GOTO 40
      H(2)=U-X
      S=U+X
      DO 20 I=3,N1
  20  H(I)=(T-S/(I-1))*H(I-1)-(U-U/(I-1))*H(I-2)
  40  CONTINUE
      RETURN
      END
