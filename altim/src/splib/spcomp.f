**SPCOMP -- Single Discrete Fourier Transform (DFT) component of Real data
*+
      COMPLEX*16 FUNCTION SPCOMP (X, N, F)
      INTEGER*4  N
      REAL*8     X(0:N-1), F
*
* This routine computes a single DFT component of a series of N equally spaced
* real data points X(0) through X(N-1).
* The Discrete Fourier Transform of X for one specific frequency F (with 1 being
* the sampling frequency) is a COMPLEX*16 value SPCOMP, whose real part is the
* cosine component and imaginary part is the sine component.
* The number of data points may be even or odd, but must be at least 2.
*
* Arguments:
*   X       (input): Real values of the N data points.
*   N       (input): Number of data points.
*   F       (input): Frequency associated with the DFT (measured in sampling
*                    frequencies).
*   SPCOMP (output): Discrete Fourier Transform of X for frequency F.
*-
* This subroutine was taken from "Signal Processing Algorithms;
* by Samuel D. Stearns and Ruth A. David, Prentice-Hall Inc., Englewood Cliffs,
* New Jersey, 1988" and was adjusted by Remko Scharroo, DUT/SOM.
*-
* 11-Mar-1991: Created, new manual, double precision implemented.
* 11-Jan-1994: Standardize.
*------------------------------------------------------------------------------
      INTEGER*4 K
      COMPLEX*16 TPJF
      TPJF=DCMPLX(0D0,-8*ATAN(1D0)*F)
      SPCOMP=X(0)
      DO K=1,N-1
         SPCOMP=SPCOMP+X(K)*EXP(K*TPJF)
      ENDDO
      END
