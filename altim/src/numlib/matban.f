C    ****************************************************************** 
C                                                      38. MATBAN
C                                                                       
C    PURPOSE                                                            
C    TO TRANSFORM A BAND MATRIX A INTO A FULL MATRIX C
C    SUITED TO BE USED IN SUBROUTINE LINBAN
C                                                                       
C    USAGE                                                              
C    CALL MATBAN(M,N,K,L,A,C)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    M    I WIDTH OF THE BAND OF A
C    N    I NUMBER OF EQUATIONS                                         
C    K    I NUMBER OF LOWER CODIAGONALS
C    L    I FIRST DIMENSION OF A AS DECLARED IN THE CALLING             
C           (SUB)PROGRAM (L.GE.N)                                       
C    A    I REAL ARRAY OF DIMENSION (L,M) SEE REMARKS.
C    C    O RESULTANT MATRIX OF DIMENSION (L,M) TO BE USED IN LINBAN
C                                                                       
C    *******************************************************************
      SUBROUTINE MATBAN(M,N,K,L,A,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(L,N),C(L,M)
      IJ=M-K-1
      IX=0
      DO 80 I=1,N
      IF (I.LE.K+1) THEN
      IJ=IJ+1
      ELSE
      IX=IX+1
      END IF
      IF (IX+IJ.GT.N) IJ=IJ-1
      DO 40 J=1,IJ
  40  C(I,J)=A(I,IX+J)
  80  CONTINUE
      RETURN
      END
