**SPCDEMN -- Complex demodulation of a non-equally spaced series of Real data
*+
      SUBROUTINE SPCDEMN (N, T, X, Y, W, F, NSMTH)
      INTEGER*4  N, NSMTH
      REAL*8     T(0:N-1), X(0:N-1), F
      COMPLEX*16 Y(0:N-1), W(0:N-1)
*  or REAL*8     Y(0:2*N-1), W(0:2*N-1)
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
*   W            : Working space.
*   F     (input): Frequency associated with the signal to be removed
*                  (measured in sampling frequencies).
*   NSMTH (input): Smoothing factor (default=1)
*-
* 25-Mar-1993: Created from SPCDEM.
* 11-Jan-1994: Standardize.
*------------------------------------------------------------------------------
      INTEGER*4  N1, K, M, M0, M1, MT
      COMPLEX*16 TPJF, SUM
      REAL*8     RUN, HRUN
*
* Set TPJF to 2*PI*J*F  (J**2=-1).
* Compute length of running average filter.
*
      TPJF=DCMPLX(0D0,-8*ATAN(1D0)*F)
      RUN=1/F
      IF (NSMTH.GT.1) RUN=NSMTH/F
      HRUN=RUN/2
      N1=N-1
*
* Multiply the data by exp(j*omega*t)
*
      DO K=0,N1
         W(K)=X(K)*EXP(T(K)*TPJF)
      ENDDO
*
* Now loop through all data.
*
      M0=0
      M1=-1
      SUM=0
      DO K=0,N1
*
* Compute running average around each data point, in three steps:
* 1. Add points to the average up to HRUN from the data point
*
         MT=M1
         DO M=MT+1,N1
            IF (T(M)-T(K).LE.HRUN) THEN
               SUM=SUM+W(M)
               M1=M
            ELSE
               GOTO 40
            ENDIF
         ENDDO
*
* 2. Drop points more than RUN before the end of the bin
*
   40    MT=M0
         DO M=MT,M1
            IF (T(M1)-T(M).GT.RUN) THEN
               SUM=SUM-W(M)
               M0=M+1
            ELSE
               GOTO 60
            ENDIF
         ENDDO
*
* 3. Include points less than RUN from the beginning of the bin
*
   60    MT=M1
         DO M=MT+1,N1
            IF (T(M)-T(M0).LE.RUN) THEN
               SUM=SUM+W(M)
               M1=M
            ELSE
               GOTO 80
            ENDIF
         ENDDO
*
* Store twice the average over the bin in Y. This is the modulated
* signal amplitude.
*
   80    CONTINUE
         Y(K)=2*SUM/(M1-M0+1)
      ENDDO
*
* The demodulated signal X(K) is computed by subtracting
* the signal determined by Y(K).
*
      DO K=0,N1
         X(K)=DBLE(X(K)-Y(K)*EXP(-T(K)*TPJF))         
      ENDDO
      END
