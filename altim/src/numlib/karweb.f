C    ****************************************************************** 
C                                                      26. KARWEB
C                                                                       
C    PURPOSE                                                            
C    TO COMPUTE THE CHARCTERISTIC EQUATION OF A ARBITRARY
C    REAL MATRIX BY MEANS OF THE WEBER-VOETER METHOD
C                                                                       
C    USAGE                                                              
C    CALL KARWEB(N,K,IA,A,V,EP)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I ORDER OF THE MATRIX A
C    K    O DEGREE OF THE COMPUTED POLYNOMIAL, IN GENERAL K=N
C           IF DURING THE PROCESS PIVOT ELEMENTS BECOME
C           LESS EP, THEN K<N
C    IA   I THE FIRST DIMENSION OF A, AS DECLARED IN THE
C           CALLING (SUB)PROGRAM. (IA.GE.N)
C    A    I THE GIVEN REAL MATRIX A(N,N), DESTROYED ON EXIT
C    V    O V(N+1) CONTAINS THE RESULT POLYNOMIAL PK(X)
C           PK(X)=V(K+1)*X**K+.......+V(1)*X**0 AND
C           V(J)=0 FOR J=K+2.......N+1
C    EP   I TOLERANCE VALUE FOR PIVOT ELEMENTS
C
C     METHODS
C     THE WEBER-VOETER METHOD, SEE
C     HOUSEHOLDER, A.S. "THE THEORY OF MATRICES IN NUMERICAL ANALYSIS"
C
C    *******************************************************************
      SUBROUTINE KARWEB(N,K,IA,A,V,EP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,N),V(*)
      DATA Z,U/0.D0,1.D0/
      DO 120 I=1,N-1
      T=Z
      M=I+1
      DO 10 L=M,N
      S=DABS(A(L,I))
      IF (S.GT.T) THEN
      T=S
      J=L
      ENDIF
   10 CONTINUE
      IF (T.LT.EP) THEN
      K=I
      GOTO 130
      ENDIF
      IF (J.NE.M) THEN
      DO 20 L=1,N
      T=A(M,L)
      A(M,L)=A(J,L)
      A(J,L)=T
      T=A(L,M)
      A(L,M)=A(L,J)
   20 A(L,J)=T
      T=A(M,M)
      A(M,M)=A(J,M)
      A(J,M)=T
      T=A(M,M)
      A(M,M)=A(M,J)
      A(M,J)=T
      ENDIF
      DO 30 L=M,N
   30 V(L+1)=A(L,I)
      DO 60 L=1,I
      T=Z
      DO 40 J=1,L-1
   40 T=T+A(L+1,J)*V(J+1)
      DO 50 J=M,N
   50 T=T+A(L+1,J)*V(J+1)
   60 V(L+1)=-T/A(L+1,L)
      DO 90 L=1,M
      IF (L.EQ.1) THEN
      T=Z
      DO 70 J=M,N
   70 T=T+A(1,J)*V(J+1)
      ELSE
      T=A(L-1,I)-V(L)
      ENDIF
      DO 80 J=L,I
   80 T=T+A(L,J)*V(J+1)
   90 A(L,M)=T
      DO 110 L=I+2,N
      T=Z
      DO 100 J=1,N
  100 T=T+A(L,J)*V(J+1)
  110 A(L,M)=T
  120 CONTINUE
      K=N
  130 CONTINUE
      DO 140 I=K+2,N+1
  140 V(I)=Z
      DO 150 I=1,K
  150 V(I)=-A(I,K)
      V(K+1)=U
      RETURN
      END
