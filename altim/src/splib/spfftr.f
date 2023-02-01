**SPFFTR -- Forward/Reverse Fast Fourier Transform (FFT) of Real data.
*+
      SUBROUTINE SPFFTR (X, N, ISIGN)
      INTEGER*4  N, ISIGN
*     REAL*8     X(0:N+1)
*
* This routine computes the FFT of a series of N [=2**k] equally spaced real
* data points X(0) through X(N-1), or it performs the reverse operation.
* The Fast Fourier Transform of X is stored back into array X. The elements
* X(0) through X(N+1) are the respective cosine and sine components of the FFT,
* starting with the zero-frequency pair [X(0),X(1)], up to the Nyquist frequency
* components [X(N),X(N+1)].
* Use ISIGN=-1 for the forward FFT; ISIGN=+1 for the reverse FFT.
* When the reverse FFT is computed, the Fourier transform is transformed back
* into N data points X scaled by a factor N.
* The number of data points must be a power of 2, i.e. 2, 4, 8, etc.
*
* Arguments:
*   N      (input): Number of data points.
*   ISIGN  (input): = -1 for forward FFT.
*       X  (input): Real values of the N data points.
*         (output): Fast Fourier Transform of X.
*   ISIGN  (input): = +1 for reverse FFT.
*       X  (input): Fast Fourier Transform.
*         (output): Real values of the N data points scaled by a factor N.
*-
* This subroutine was taken from "Signal Processing Algorithms;
* by Samuel D. Stearns and Ruth A. David, Prentice-Hall Inc., Englewood Cliffs,
* New Jersey, 1988" and was adjusted by Remko Scharroo, DUT/SOM.
*-
* 11-Mar-1991: Created, new manual, double precision implemented, forward and
*              inverse FFT combined.
* 11-Jan-1994: Standardize.
*------------------------------------------------------------------------------
      INTEGER M
      COMPLEX*16 X(0:N/2),U,D,TMP
      REAL*8 TPN
      TPN=8*ATAN(1D0)/N
      IF (ISIGN.LT.0) THEN
         CALL SPFFTC(X,N/2,-1)
         X(N/2)=X(0)
         DO M=0,N/4
            U=DCMPLX(SIN(M*TPN),COS(M*TPN))
            D=DCONJG(X(N/2-M))
            TMP =((1+U)*X(M)+(1-U)*D)/2
            X(M)=((1-U)*X(M)+(1+U)*D)/2
            X(N/2-M)=DCONJG(TMP)
         ENDDO
      ELSE
         DO M=0,N/4
            U=DCMPLX(SIN(M*TPN),-COS(M*TPN))
            D=DCONJG(X(N/2-M))
            TMP =(1+U)*X(M)+(1-U)*D
            X(M)=(1-U)*X(M)+(1+U)*D
            X(N/2-M)=DCONJG(TMP)
         ENDDO
         CALL SPFFTC(X,N/2,+1)
      ENDIF
      END
