**SPFAMP -- Forward/Reverse change between Fourier Transform and Series.
*+
      SUBROUTINE SPFAMP (X, N, ISIGN)
      INTEGER*4  N, ISIGN
      COMPLEX*16 X(0:N/2)
*  or REAL*8     X(0:N+1)
*
* This routine computes the cosine and sine amplitudes for each frequency
* component of a Fourier series out of the Fourier Transform, or vice versa.
* The array X is either COMPLEX*16 or REAL*8. When interpreted as REAL*8
* X(0) through X(N+1) are the respective cosine and sine components of the
* Fourier Transform, or the cosine and sine amplitudes of the Fourier Series.
*
* Arguments:
*   N      (input): Number of data points.
*   ISIGN  (input): = -1 for forward transformation.
*       X  (input): Fourier Transform components.
*         (output): Fourier Series amplitudes.
*   ISIGN  (input): = +1 for reverse transformation.
*       X  (input): Fourier Series amplitudes.
*         (output): Fourier Transform components.
*-
* 11-Mar-1991: Created.
* 11-Jan-1994: Standardize.
*------------------------------------------------------------------------------
      INTEGER M
      REAL*8 FAC
      FAC=(N/1D0)**ISIGN
      X(0)=FAC*DCONJG(X(0))
      FAC=(N/2D0)**ISIGN
      DO M=1,N/2
         X(M)=FAC*DCONJG(X(M))
      ENDDO
      END
