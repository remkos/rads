C    ****************************************************************** 
C                                                      40. MATBMV
C                                                                       
C    PURPOSE                                                            
C    TO MULTIPLY A VECTOR X BY A BAND MATRIX A (AX=B), WHERE
C    MATRIX A IS SUITED TO BE USED IN SUBROUTINE LINBAN
C                                                                       
C    USAGE                                                              
C    CALL MATBMV(M,N,K,L,A,X,B)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    M    I WIDTH OF THE BAND OF A
C    N    I NUMBER OF EQUATIONS                                         
C    K    I NUMBER OF LOWER CODIAGONALS
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING             
C           (SUB)PROGRAM (L.GE.N)                                       
C    A    I REAL ARRAY OF DIMENSION (L,M) SEE REMARK 1 OF LINBAN
C    X    I THE GIVEN VECTOR
C    B    O RESULTANT VECTOR
C                                                                       
C    *******************************************************************
      SUBROUTINE MATBMV(M,N,K,L,A,X,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(L,M),X(N),B(N)
      IJ=M-K-1
      IX=0
      DO 80 I=1,N
      IF (I.LE.K+1) THEN
      IJ=IJ+1
      ELSE
      IX=IX+1
      END IF
      IF (IX+IJ.GT.N) IJ=IJ-1
      T=0.0D0
      DO 40 J=1,IJ
  40  T=T+A(I,J)*X(IX+J)
  80  B(I)=T
      RETURN
      END
