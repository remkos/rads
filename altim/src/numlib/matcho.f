C    ****************************************************************** 
C                                                      41. MATCHO
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE CHOLESKI REDUCTION A=LDU OF A REAL
C    MATRIX A, WITHOUT CALCULATING SQUARE ROOTS,
C    WHERE A IS A SYMMETRIC MATRIX IN VECTOR FORM
C                                                                       
C    USAGE                                                              
C    CALL MATCHO(A,H,V,KP,N,D,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I REAL VECTOR CONTAINING THE COLUMS OF THE UPPER PART OF
C           THE SYMMETRIC MATRIX A. THE DIMENSION IS AT LEAST N*(N+1)/2
C           A CAN BE FORMED BY SUBROUTINE MATSYD
C         O IF IF=0, A CONTAINS THE DIAGONAL MATRIX D AND THE
C           UPPER DIAGONAL MATRIX U OF THE REDUCTION A=LDU
C           IF IF=1, D AND U ARE PERMUTATED BY SEARCHING THE GREATEST
C           PIVOT ON THE MAIN DIAGONAL, TO  BE USED IN LINCHO AND
C           INVCHO
C    H      INTEGER VECTOR WITH H(I)=I*(I-1)/2
C           H CAN BE FILLED WITH MATSY1
C    V      WORKING SPACE OF DIMENSION AT LEAST N
C    KP   O PERMUTATION VECTOR, TO BE USED IN LINCHO AND INVCHO
C    N    I DIMENSION OF SYMMETRIX MATRIX A IS (N,N)
C    D    O DETERMINANT OF A
C    IF   I INTEGER SWITCH
C           IF=0 THE NORMAL REDUCTION A=LDU IS COMPUTED
C           IF=1 PIVOTTING WILL TAKE PLACE TO BE USED IN LINCHO
C           AND INVCHO
C         O INTEGER ERRORCODE
C           IF=0 NORMAL EXIT
C           IF=1 PIVOT = 0
C
C     REMARKS
C     MATCHO USES SUBROUTINE MATWIS TO CHANGE TWO ROWS
C     AND THE CORRESPONDING COLUMS IN MATRIX A
C
* Debugging by Remko Scharroo, Delft, December 17, 1990
* 16-Feb-1998 - Replaced some J*(J-1)/2 by H(J)
C    *******************************************************************
      SUBROUTINE MATCHO(A,H,V,KP,N,D,IF)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER H(N)
      DIMENSION A(*),V(N),KP(N)
      DATA Z,U/0.D0,1.0D0/
      D=U
      DO 20 K=1,N
         KP(K)=K
         IF (IF.EQ.1) THEN
            T=Z
            DO 4 J=K,N
               S=DABS(A(H(J)+J))
               IF (S.GT.T) THEN
                  T=S
                  I=J
               ENDIF
   4        CONTINUE
            KP(K)=I
            IF (K.NE.I) CALL MATWIS(A,H,N,K,I)
         ENDIF
         T=A(H(K)+K)
         D=D*T
         IF (T.EQ.Z) GOTO 30
         DO 10 J=K+1,N
            M=K+H(J)
            V(J)=A(M)
  10        A(M)=V(J)/T
         DO 20 J=K+1,N
            M=H(J)
            T=A(K+M)
            DO 20 I=K+1,J
  20           A(I+M)=A(I+M)-V(I)*T
      IF=0
      RETURN
  30  IF=1
      END

      SUBROUTINE MATWIS(A,H,N,K,I)
      REAL*8 A(*),T
      INTEGER H(N)
      KK=H(K)
      II=H(I)
      DO 10 J=1,K-1
         T=A(J+KK)
         A(J+KK)=A(J+II)
  10     A(J+II)=T
      DO 20 J=K+1,I-1
         JJ=H(J)
         T=A(K+JJ)
         A(K+JJ)=A(J+II)
  20     A(J+II)=T
      DO 30 J=I+1,N
         JJ=H(J)
         T=A(K+JJ)
         A(K+JJ)=A(I+JJ)
  30     A(I+JJ)=T
      T=A(K+KK)
      A(K+KK)=A(I+II)
      A(I+II)=T
      END
