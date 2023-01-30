**SPDFTC -- Forward/Reverse Discrete Fourier Transform (DFT) of Complex data.
*+
      SUBROUTINE SPDFTC (X, Y, N, ISIGN)
      INTEGER    N, ISIGN
      COMPLEX*16 X(0:N-1), Y(0:N-1)
*
* This routine computes the DFT of a series of N equally spaced complex data
* points X(0) through X(N-1), or the performs the reverse operation.
* The Discrete Fourier Transform of the COMPLEX*16 array X is stored in the
* COMPLEX*16 array Y, which may not be the same as array X.
* The elements Y(0) through Y(N-1) denote the components of the DFT, from
* frequency = 0 up to (N-1)/N.
* Use ISIGN=-1 for the forward DFT, ISIGN=+1 for the reverse DFT.
* When the reverse DFT is computed, Y is transformed back into N data points X
* scaled by a factor N.
* The number of data points may be even or odd, but must be at least 2.
*
* Arguments:
*   N      (input): Number of data points.
*   ISIGN  (input): = -1 for forward DFT.
*       X  (input): Complex values of the N data points.
*       Y (output): Discrete Fourier Transform of X.
*   ISIGN  (input): = +1 for reverse DFT.
*       X (output): Complex values of the N data points scaled by a factor N.
*       Y  (input): Discrete Fourier Transform of X.
*-
* This subroutine was taken from "Signal Processing Algorithms;
* by Samuel D. Stearns and Ruth A. David, Prentice-Hall Inc., Englewood Cliffs,
* New Jersey, 1988" and was adjusted by Remko Scharroo, DUT/SOM.
*-
* 11-Mar-1991: Created, new manual, double precision implemented.
*  2-Jul-1993: Standardize manual.
*------------------------------------------------------------------------------
      INTEGER*4 K,M
      COMPLEX*16 TPJN
      TPJN=DCMPLX(0D0,8*ISIGN*ATAN(1D0)/N)
      DO M=0,N-1
         Y(M)=X(0)
         DO K=1,N-1
            Y(M)=Y(M)+X(K)*EXP(TPJN*K*M)
         ENDDO
      ENDDO
      END
