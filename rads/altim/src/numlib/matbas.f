C    ****************************************************************** 
C                                                      39. MATBAS
C                                                                       
C    PURPOSE                                                            
C    TO GENERATE AN INDEPENDENT BASE OF ROWS IN A MATRIX A
C                                                                       
C    USAGE                                                              
C    CALL MATBAS(A,JR,B,P,Q,IA,M,IB,N,EP,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I M BY N COEFFICIENT MATRIX, UNCHANGED ON EXIT
C    JR   O THE RESULT NUMBERS OF THE ROWS IN A,
C           WHICH FORM AN INDEPENDENT BASE
C    B    O N BY N AUXILIARY WORKING ARRAY, ON EXIT B
C           CONTAINS THE RESULT BASE ROWS IN ORTHOGONAL FORM
C    P    O AN AUXILIARY WORKING VECTOR OF N ELEMENTS, ON EXIT
C           P CONTAINS THE SQUARES OF THE LENGTHS
C           OF THE THE ROWS IN B
C    Q      AN AUXILIARY WORKING VECTOR OF N ELEMENTS
C    IA   I THE FIRST DIMENSION OF A, AS DECLARED IN THE
C           CALLING (SUB)PROGRAM. (IA.GE.M)
C    M    I FIRST DIMENSION OF A
C    IB   I THE FIRST DIMENSION OF B, AS DECLARED IN THE
C           CALLING (SUB)PROGRAM. (IB.GE.N)
C    N    I FIRST DIMENSION OF B
C    EP   I UPPER BOUND OF THE TOLERANCE FOR THE LENGTHS
C           OF THE BASE ROWS
C    IF   O INTEGER ERRORCODE
C           IF=0 NORMAL EXIT, AN INDEPENDENT BASE OF N ROWS OF
C                A IS FOUND AND THE ROW NUMBERS ARE GIVEN IN JR
C           IF=1 NO N INDEPENDENT ROWS ARE FOUND USING THE GIVEN
C                VALUE OF EP. IN JR THE ROWNUMBERS ARE GIVEN OF
C                THE ROWS IN A, WHICH FORM AN INDEPENDENT BASE
C                OF RANK LESS THAN N
C
C     METHODS
C     THE GRAMM-SCHMIDT ORTHOGONALISATION PROCESS
C
C    *******************************************************************
      SUBROUTINE MATBAS(A,JR,B,P,Q,IA,M,IB,N,EP,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,N),JR(N),B(IB,N),P(N),Q(N)
      DATA Z/0.D0/
      IF=0
      DO 20 I=1,M
      T=Z
      II=I
      DO 10 J=1,N
      S=A(I,J)
      B(1,J)=S
   10 T=T+S*S
      IF (T.GT.EP) GOTO 30
   20 CONTINUE
      GOTO 100
   30 P(1)=T
      JR(1)=II
      KK=1
      DO 90 I=II+1,M
      T=Z
      DO 40 J=1,N
   40 T=T+A(I,J)*A(I,J)
      R=T
      DO 60 K=1,KK
      S=Z
      DO 50 J=1,N
   50 S=S+A(I,J)*B(K,J)
      Q(K)=S/P(K)
      R=R-S*Q(K)
   60 CONTINUE
      IF (R.LT.T*EP) GOTO 90
      KK=KK+1
      DO 80 J=1,N
      T=A(I,J)
      DO 70 K=1,KK-1
   70 T=T-Q(K)*B(K,J)
   80 B(KK,J)=T
      P(KK)=R
      JR(KK)=I
      IF (KK.EQ.N) GOTO 110
   90 CONTINUE
  100 IF=1
  110 CONTINUE
      RETURN
      END
