C    ****************************************************************** 
C                                                      33. LINCHO
C                                                                       
C    PURPOSE                                                            
C    TO SOLVE A SYSTEM OF LINEAR EQUATIONS AX=B BY CHOLESKI'S METHOD
C    WHERE A IS A SYMMETRIC MATRIX IN VECTOR FORM
C                                                                       
C    USAGE                                                              
C    CALL LINCHO(A,H,X,V,KP,N,D,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX A. THE DIMENSION IS AT LEAST N*(N+1)/2
C           A CAN BE FORMED BY SUBROUTINE MATSYD
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH MATSY1
C    X    I RIGHT HAND SIDE B
C         O SOLUTION VECTOR X
C    V      WORKING SPACE OF DIMENSION AT LEAST N
C    KP     WORKING SPACE OF DIMENSION AT LEAST N
C    N    I DIMENSION OF SYMMETRIX MATRIX A IS (N,N)
C    D    O DETERMINANT OF A
C    IF   O INTEGER ERRORCODE
C           IF=0 NORMAL EXIT
C           IF=1 PIVOT = 0
C
C     REMARKS
C     LINCHO USES SUBROUTINE MATCHO and LINCH1
C
C    *******************************************************************
      SUBROUTINE LINCHO(A,H,X,V,KP,N,D,IF)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER H(N)
      DIMENSION A(*),X(N),V(N),KP(N)
      IF=1
      CALL MATCHO(A,H,V,KP,N,D,IF)
      IF (IF.EQ.1) GOTO 50
      CALL LINCH1(A,H,X,KP,N)
   50 END
