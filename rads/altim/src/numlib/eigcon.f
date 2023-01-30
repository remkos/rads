C    ****************************************************************** 
C                                                       8. EIGCON
C                                                                       
C    PURPOSE                                                            
C    TO CHECK THE ACCURACY OF THE EIGENVALUES AND EIGENVECTORS
C    PROVIDED BY SUBROUTINE EIGHES BY CALCULATING A*X-L*X
C                                                                       
C    USAGE                                                              
C    CALL EIGCON(N,K,IA,IE,A,EV,EW,U)
C                                                                       
C    DESCRIPTION OF PARAMETERS                                          
C    N    I ORDER OF THE MATRIX A
C    K    I NUMBER OF EIGENVALUES AND EIGENVECTORS TO BE CHECKED
C    IA   I FIRST DIMENSION OF A AS DECLARED IN THE
C           DIMENSION STATEMENT OF THE CALLING PROGRAM
C    IE   I FIRST DIMENSION OF EV AND EW AS DECLARED IN THE
C           DIMENSION STATEMENT OF THE CALLING PROGRAM
C    A    I ARRAY A(IA,N) CONTAINS THE MATRIX A.
C    EW   I ARRAY EW(IE,IE,2) CONTAINS THE EIGENVECTORS X
C    EV   I ARRAY EV(IE,2) CONTAINS THE EIGENVALUES L
C    U    0 RESULTANT VECTOR WITH THE MAXIMUM ERRORS, SEE REMARK 1
C
C    REMARKS
C    1. VECTOR U IS FILLED WITH THE MAXIMUM OF THE 2*N VALUES
C       OF THE COMPLEX VECTORS A*X-L*X FOR EACH EIGENVECTOR X
C       AND EIGENVALUE L.
C
C    *******************************************************************
      SUBROUTINE EIGCON(N,K,IA,IE,A,EV,EW,U)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,N),EV(IE,IE,2),EW(IE,2),U(N)
      ZZ=0.D0
      DO 30 J=1,K
      X=ZZ
      Y=ZZ
      DO 20 I=1,N
      S=ZZ
      T=ZZ
      DO 10 M=1,N
      T=T+A(I,M)*EV(M,J,1)
   10 S=S+A(I,M)*EV(M,J,2)
      T=T-EW(J,1)*EV(I,J,1)+EW(J,2)*EV(I,J,2)
      S=S-EW(J,1)*EV(I,J,2)-EW(J,2)*EV(I,J,1)
      IF (X.LT.DABS(T)) X=DABS(T)
      IF (Y.LT.DABS(S)) Y=DABS(S)
   20 CONTINUE
      U(J)=X
      IF (X.LT.Y) U(J)=Y
   30 CONTINUE
      RETURN
      END
