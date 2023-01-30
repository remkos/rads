C    ****************************************************************** 
C                                                      21. INVCHO
C                                                                       
C    PURPOSE                                                            
C    TO INVERT A MATRIX A BY CHOLESKI'S REDUCTION
C    WHERE A IS A SYMMETRIC MATRIX IN VECTOR FORM
C                                                                       
C    USAGE                                                              
C    CALL INVCHO(A,H,V,KP,N,D,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX A. THE DIMENSION IS AT LEAST N*(N+1)/2
C           A CAN BE FORMED BY SUBROUTINE MATSYD
C         O INVERTED MATRIX A, STORED AS DESCRIBED ABOVE
C    H    I INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH MATSY1
C    V      WORKING SPACE OF DIMENSION AT LEAST N
C    KP     WORKING SPACE OF DIMENSION AT LEAST N
C    N    I DIMENSION OF SYMMETRIX MATRIX A IS (N,N)
C    D    O DETERMINANT OF A
C    IF   O INTEGER ERRORCODE
C           IF=0 NORMAL EXIT
C           IF=1 PIVOT = 0
C
C     REMARKS
C     INVCHO USES THE SUBROUTINES MATCHO AND MATWIS
C
C Debugged by Remko Scharroo, Delft, December 14, 1990.
C    *******************************************************************
      SUBROUTINE INVCHO(A,H,V,KP,N,D,IF)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER H(N)
      DIMENSION A(*),V(N),KP(N)
      IF=1
      CALL MATCHO(A,H,V,KP,N,D,IF)
      IF (IF.EQ.1) GOTO 100
      DO 30 J=2,N
         J1=J-1
         JJ=H(J)
         DO 10 K=2,J1
   10       V(K)=A(K+JJ)
         DO 30 I=J1,1,-1
            T=A(I+JJ)
            DO 20 K=I+1,J1
   20          T=T+V(K)*A(I+H(K))
   30       A(I+JJ)=-T
      DO 50 I=1,N
         II=I*(I+1)/2
   50    A(II)=1/A(II)
      DO 70 K=2,N
         KK=H(K)
         DO 60 I=K-1,1,-1
            V(I)=A(I+KK)
   60       A(I+KK)=V(I)*A(K+KK)
         DO 70 J=1,K-1
            JJ=H(J)
            DO 70 I=1,J
   70          A(I+JJ)=A(I+JJ)+V(J)*A(I+KK)
      DO 90 K=N,1,-1
         J=KP(K)
   90    IF (K.NE.J) CALL MATWIS(A,H,N,K,J)
  100 END
