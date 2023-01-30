**SPFFTC -- Forward/Reverse Fast Fourier Transform (FFT) of Complex data.
*+
      SUBROUTINE SPFFTC (X, N, ISIGN)
      INTEGER*4  N, ISIGN
      COMPLEX*16 X(0:N-1)
*
* This routine computes the FFT of a series of N equally spaced complex data
* points X(0) through X(N-1), or the performs the reverse operation.
* The Fast Fourier Transform of the COMPLEX*16 array X is stored back into the
* array X itself. The elements X(0) through X(N-1) denote the components of the
* FFT, from frequency = 0 up to (N-1)/N.
* Use ISIGN=-1 for the forward FFT, ISIGN=+1 for the reverse FFT.
* When the reverse FFT is computed, the Fourier transform of X is transformed
* back into N data points X scaled by a factor N.
* The number of data points must be a power of 2, i.e. 2, 4, 8, etc.
*
* Arguments:
*   N      (input): Number of data points.
*   ISIGN  (input): = -1 for forward FFT.
*       X  (input): Complex values of the N data points.
*         (output): Fast Fourier Transform of X.
*   ISIGN  (input): = +1 for reverse FFT.
*       X  (input): Fast Fourier Transform.
*         (output): Complex values of the N data points scaled by a factor N.
*-
* This subroutine was taken from "Signal Processing Algorithms;
* by Samuel D. Stearns and Ruth A. David, Prentice-Hall Inc., Englewood Cliffs,
* New Jersey, 1988" and was adjusted by Remko Scharroo, DUT/SOM.
*-
* 11-Mar-1991: Created, new manual, double precision implemented.
* 11-Jan-1994: Standardize.
*------------------------------------------------------------------------------
      INTEGER I,L,M,MR
      COMPLEX*16 T,PJ
      PJ=DCMPLX(0D0,4*ISIGN*ATAN(1D0))
      MR=0
      DO M=1,N-1
         L=N
    1    L=L/2
         IF (MR+L.GE.N) GOTO 1
         MR=MOD(MR,L)+L
         IF (MR.LE.M) GOTO 2
         T=X(M)
         X(M)=X(MR)
         X(MR)=T
    2    continue
      ENDDO
      L=1
    3 IF (L.GE.N) RETURN
      DO M=0,L-1
         DO I=M,N-1,2*L
            T=X(I+L)*EXP(M*PJ/L)
            X(I+L)=X(I)-T
            X(I)=X(I)+T
         ENDDO
      ENDDO
      L=2*L
      GOTO 3
      END
