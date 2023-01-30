C    ****************************************************************** 
C                                                      27. KLEFUN
C                                                                       
C    PURPOSE                                                            
C    TO ESTIMATE THE VALUES OF THE PARAMETERS IN A GIVEN FUNCTION
C    FUN WHICH IS APPROXIMATING A SET OF (X,Y) POINTS, USING
C    THE LEAST SQUARES METHOD
C                                                                       
C    USAGE                                                              
C    CALL KLEFUN(N,M,L,EP,R,X,Y,C,FUN,DIF,V,A,DC,IF)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I THE NUMBER OF POINTS
C    M    I THE NUMBER OF PARAMETERS TO BE ESTIMATED
C    L    I FIRST DIMENSION OF A AND DC AS DECLARED
C           IN THE CALLING PROGRAM
C    EP   I RELATIVE ACCURACY PRESCIBED FOR THESE PARAMETERS
C    R    O SUM (J:1,N) (Y(J)-FUN(X(J)))**2
C    X    I X(N) CONTAINS THE ABSCISSA OF THE GIVEN POINTS
C    Y    I Y(N) CONTAINS THE ORDINATES OF THE GIVEN POINTS
C    C    I C(M) MUST CONTAIN A REASONABLE ESTIMATION
C           OF THE WANTED PARAMETERS
C         O THE COMPUTED VALUES OF THE PARAMETERS
C    FUN  I THE GIVEN FUNCTION, SEE REMARKS
C    DIF  I A FUNCTION, YIELDING THE M PARTIAL DERIVATIVES
C           DF/DC(I), SEE REMARKS
C    V      V(M) WORKING AREA OF DIMENSION OF AT LEAST M
C    A      A(L,M) WORKING AREA, WITH L AT LEAST M
C    DC     DC(L,N) WORKING AREA, WITH L AT LEAST M
C    IF   O INTEGER ERROR CODE
C                IF=0 NORMAL EXIT
C                IF=1 PIVOT LESS EP IN LINGAU
C                IF=2 NUMBER OF ITERATION STEPS EXCEEDS 15
C
C    REMARKS
C    KLEFUN USES THE FOLLOWING SUBROUTINES
C    1. SUBROUTINE LINGAU
C    2. DOUBLE PRECISION FUNCTION FUN(X,C)
C       STATEMENTS TO COMPUTE THE FUNCTION VALUE IN POINT X
C    3  SUBROUTINE DIF(J,X,C,DC)
C       STATEMENTS TO COMPUTE THE M PARTIAL DERIVATIVES DF/DC(J)
C       IN A POINT X, WHICH IS ESTABLISHED IN KLEFUN BY
C           DO 100 J=1,N
C       100 CALL DIF(J,X(J),C,DC)
C
C    REFERENCES
C    RC-TWA-75003 NUMLIBDA METHODEBESCHRIJVINGEN
C
C    *******************************************************************
      SUBROUTINE KLEFUN(N,M,L,EP,R,X,Y,C,FUN,DIF,V,A,DC,IF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),C(M),V(M),A(L,M),DC(L,N)
      EXTERNAL FUN,DIF
      DATA Z/0.D0/
      IF=0
      NT=0
   10 IF (NT.GT.15) THEN
      IF=2
      GOTO 200
      ENDIF
      DO 20 J=1,N
   20 CALL DIF(J,X(J),C,DC)
      DO 50 I=1,M
      V(I)=Z
      DO 40 J=1,N
   40 V(I)=V(I)+(Y(J)-FUN(X(J),C))*DC(I,J)
   50 CONTINUE
      DO 80 I=1,M
      DO 70 K=1,M
      A(I,K)=Z
      DO 60 J=1,N
   60 A(I,K)=A(I,K)+DC(I,J)*DC(K,J)
   70 CONTINUE
   80 CONTINUE
      CALL LINGAU(M,L,A,V,D,EP,IF)
      IF (IF.EQ.1) GOTO 200
      DO 90 I=1,M
   90 C(I)=C(I)+V(I)
      DO 100 I=1,M
      IF (DABS(V(I)).GT.EP*DABS(C(I))) THEN
      NT=NT+1
      GOTO 10
      ENDIF
  100 CONTINUE
      R=Z
      DO 110 I=1,N
  110 R=R+(Y(I)-FUN(X(I),C))**2
  200 CONTINUE
      RETURN
      END
