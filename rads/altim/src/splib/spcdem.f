**SPCDEM -- Complex demodulation of a series of Real data
*+
      SUBROUTINE SPCDEM (N, X, Y, F, NSMTH)
      INTEGER    N, NSMTH
      REAL*8     X(0:N-1), F
      COMPLEX*16 Y(0:N-1)
*  or REAL*8     Y(0:2*N-1)
*
* This routine performs a complex demodulation of a single frequency
* component in a series of N equally spaced real data points X(0) through
* X(N-1). This means that the varying amplitude and phase of the
* signal component with frequency F is determined and subtracted from the
* data. On return, Y(0) through Y(N-1) will contain the complex amplitude
* of the signal with frequency F at data points X(0) through X(N-1).
* When interpreted as REAL*8, Y(0) and Y(1) are the cosine and sine
* amplitudes in X(0), etc.
* X(0) through X(N-1) will then be replaced by the demodulated signal,
* i.e. real values with the modulated signal removed.
* The smoothing factor NSMTH determines how many cycles are used for the
* signal demodulation. Optimally this should be 1, but higher values
* are allowed to increase the smoothing.
*
* Arguments:
*   N     (input): Number of data points.
*   X     (input): Real values of the N data points.
*        (output): Demodulated signal in each of the N data points.
*   Y    (output): Modulation (complex amplitude) of the signal component with
*                  frequency F.
*   F     (input): Frequency associated with the signal to be removed
*                  (measured in sampling frequencies).
*   NSMTH (input): Smoothing factor (default=1)
*-
* 11-Mar-1993: Created.
* 11-Jan-1994: Standardize.
*------------------------------------------------------------------------------
      INTEGER*4  LRUN, LRUN1, LRUN2, K, M
      COMPLEX*16 TPJF, SUM, Y0
      REAL*8     FAC
*
* Set TPJF to 2*PI*J*F  (J**2=-1).
* Compute length of running average filter.
*
      TPJF=DCMPLX(0D0,-8*ATAN(1D0)*F)
      LRUN=NINT(1/F)
      IF (NSMTH.GT.1) LRUN=NINT(NSMTH/F)
      LRUN=MAX(2,LRUN)
      LRUN1=LRUN-1
      LRUN2=LRUN/2
      FAC=2D0/LRUN
*
* Multiply the data by exp(j*omega*t)
*
      DO K=0,N-1
         Y(K)=X(K)*EXP(K*TPJF)
      ENDDO
*
* Compute the sum over the first running average bin and store
* twice the average value in Y(0).
*
      SUM=Y(0)
      DO K=1,LRUN1
         SUM=SUM+Y(K)         
      ENDDO
      Y0=Y(0)
      Y(0)=FAC*SUM
*
* Do the same for the next running average bins by subtracting the previous
* modulated data point and adding the next.
*
      DO K=1,N-LRUN
         SUM=SUM-Y0+Y(K+LRUN1)
         Y0=Y(K)
         Y(K)=FAC*SUM
      ENDDO
*
* Y(K) contains now of the modulated signal over a bin extending
* from K to K+LRUN-1. Y(K) should be shifted forward by LRUN/2.
* And the demodulated signal X(K) is computed by subtracting
* the signal determined by Y(K).
*
      DO K=N-1,0,-1
         M=MAX(0,MIN(K-LRUN2,N-LRUN))
         Y(K)=Y(M)
         X(K)=DBLE(X(K)-Y(K)*EXP(-K*TPJF))         
      ENDDO
      END
