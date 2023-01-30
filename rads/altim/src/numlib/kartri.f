C    ****************************************************************** 
C                                                      25. KARTRI
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE CHARCTERISTIC EQUATION OF A ARBITRARY
C    REAL TRIDIAGONAL MATRIX
C                                                                       
C    USAGE                                                              
C    CALL KARTRI(N,A,C,U,L)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I ORDER OF THE MATRIX A
C    A    O A(N+1) CONTAINS THE RESULT POLYNOMIAL PN(X)
C           PN(X)=A(N+1)*X**N+.......+A(1)*X**0 AND
C    B      AUXILIARY WORKING VECTOR OF DIMENSION AT LEAST N+1
C    C    I C(N) CONTAINS THE MAINDIAGONAL ELEMENTS
C    U    I U(N-1) CONTAINS THE UPPER DIAGONAL ELEMENTS
C    L    I L(N-1) CONTAINS THE LOWER DIAGONAL ELEMENTS
C
C    REMARKS
C    IF THE GIVEN MATRIX IS SYMMETRIC, THEN U=L
C
C    *******************************************************************
      SUBROUTINE KARTRI(N,A,B,C,U,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION L
      DIMENSION A(*),B(*),C(*),U(*),L(*)
      DATA Z/0.D0/
      DO 10 I=2,N+1
      A(I)=Z
   10 B(I)=Z
      B(1)=1.D0
      A(2)=1.D0
      A(1)=-C(1)
      DO 30 I=2,N
      A(I+1)=1.D0
      B(I)=1.D0
      T=L(I-1)*U(I-1)
      A(I)=A(I-1)-C(I)
      DO 20 J=I-2,1,-1
      S=A(J+1)
      A(J+1)=A(J)-C(I)*S-T*B(J+1)
   20 B(J+1)=S
      S=A(1)
      A(1)=-C(I)*S-T*B(1)
   30 B(1)=S
      RETURN
      END
