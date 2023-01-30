C    ****************************************************************** 
C                                                      57. POLDEL
C                                                                       
C    PURPOSE                                                            
C    TO DEVIDE A POLYNOMIAL A BY ANOTHER POLYNOMIAL B
C                                                                       
C    USAGE                                                              
C    CALL POLDEL(A,B,N,M)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    A    I ARRAY A(N) CONTAINS THE COEFFICIENTS OF THE
C           POLYNOMIAL A: SUM A(J)*X**(J-1) J=1(1)N
C         O ON EXIT A CONTAINS THE COEFFICIENTS OF THE
C           RESULTANT POLYNOMIAL AND OF THE REMAINDER
C    B    I ARRAY B(M) CONTAINS THE COEFFICIENTS OF THE
C           POLYNOMIAL B: SUM B(J)*X**(J-1) J=1(1)M
C    N    I DEGREE OF POLYNOMIAL A IS (N-1)
C    M    I DEGREE OF POLYNOMIAL B IS (M-1)
C
C    REMARKS
C    1. DEGREE OF POLYNOMIALIS (M-N), SO M MUST BE LESS OR EQUAL N.
C    2. POLYNOMIAL A IS DIVIDED BY POLYNOMIAL B SUCH THAT
C       A = Q*B+R, WHERE
C       Q = SUM A(M+J)*X**(J-1)  J=1 (1) (N-M)
C       R = SUM A(J)*X**(J-1)    J=1 (1) (M-1)
C    3. B(M) IS NOT EQUAL 0
C
C    *******************************************************************
      SUBROUTINE POLDEL(A,B,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(M)
      S=B(M)
      DO 10 I=1,M-1
   10 B(I)=B(I)/S
      B(M)=1.D0
      DO 30 I=N,M,-1
      K=I-M
      DO 20 J=1,M-1
   20 A(J+K)=A(J+K)-A(I)*B(J)
   30 CONTINUE
      DO 40 I=M,N
   40 A(I)=A(I)/S
      RETURN
      END
