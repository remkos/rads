**SPECTR -- Determine spectral lines in a non-equally spaced time series.
*+
      SUBROUTINE SPECTR (N, T, X, M, F, A, B, PROB, SIGMA, L, FL, PL,
     .                   OFAC, HIFAC, FALSE, RESID, VERBOSE)
      INTEGER*4 N, M, L
      REAL*8    T(N), X(N), F(M), A(M), B(M), PROB(M), SIGMA(0:M),
     .          FL(L), PL(L), OFAC, HIFAC, FALSE, RESID
      LOGICAL   VERBOSE
*
* Compute the frequency, sine- and cosine-amplitude of the major spectral
* lines in a time series X(t). Lines that do not reach significantly above
* the continuum of the background noise will not appear in the spectrum.
* To separate the spectral lines from the continuum, the user has to specify
* a so-called `false alarm probability', i.e. the probability that lines
* from the continuum emerge in the selected spectrum. If the probability is
* set too low, physically significant peaks may be missed; if set too high,
* peaks due to noise will be mislabeled as signal peaks. Values between 0.05
* and 0.01 are suggested for this `false alarm probability', corresponding to
* one chance in 20 or one chance in 100 of mistaking a noise peak for signal,
* provided the noise is Gaussian.
* The spectral lines are isolated one by one using the Lomb-Scargle normalized
* periodogram. In every iteration the frequency of the largest peak in the
* periodogram is determined and tested against the false alarm probability. If
* the peak is significant, the amplitude and phase of the sinusoid at the peak
* frequency is determined by least squares. Afterwards the peak sinusoid is
* subtracted from the series point by point. This procedure is repeated until
* no further peaks meet the false alarm probability. Thus we have
*
*           m
*   X(t) = sum [ A  cos (2 pi f t) + B  sin (2 pi f t) ] + NOISE
*           j     j            j      j            j
*
* The frequencies that may be extracted range from zero up to HIFAC times
* the average Nyquist frequency, being 0.5/(average time interval). The
* resolution of the frequencies to be extracted can be increased by increasing
* the oversampling factor OFAC; OFAC=1 refers to a resolution equal to the Fejer
* frequency, i.e. 1/(total time interval).
*
* The iteration performed to select peaks in the periodogram is terminated
* when either:
* 1. A peak does not reach the `false alarm probability' as explained above.
* 2. The residual noise after removing the selected signals has a std.dev.
*    which is less than RESID times the std.dev. of the a priori signal.
*    (RESID has to be specified by the user.)
* 3. The maximum number of frequencies to be isolated (M) is reached.
*
* Arguments:
*  N       (input) : Number of points in the time series.
*  T       (input) : The independent variable (time).
*  X       (input) : Values in each point.
*         (output) : Residual noise after removal of all significant signals.
*  M       (input) : Maximum number of frequencies to be isolated.
*         (output) : Number of isolated frequencies that meet the false alarm
*                    probability.
*  F      (output) : Frequencies of the spectral lines in 1/(same units as T).
*  A      (output) : Cosine amplitude for each spectral line in same units as X.
*  B      (output) : Sine amplitude for each spectral line in same units as X.
*  PROB   (output) : False alarm probability of each spectral line.
*  SIGMA  (output) : Standard deviation of the signal remaining after removal
*                    of each spectral line, in the same units as X.
*                    At return SIGMA(0) will contain the A PRIORI sigma.
*  L       (input) : Dimension of the working spaces FL and PL.
*         (output) : Number of samples in the last computed periodogram.
*  FL     (output) : Frequency of each sample in this periodogram.
*  PL     (output) : Spectral density in this periodogram (dB).
*  OFAC    (input) : Oversampling factor (typical value: 4 or larger).
*  HIFAC   (input) : Maximum frequency measured in Nyquist frequencies.
*  FALSE   (input) : Limit of false alarm probability (see above).
*  RESID   (input) : Maximum fraction of noise remaining (see above).
*  VERBOSE (input) : If .TRUE., type out period, amplitude and phase (deg)
*                    during iteration.
*-
* 12-Mar-1991: Created by Remko Scharroo, DUT/SOM (c)
* 10-Dec-1992: RESID and VERBOSE added
* 14-Dec-1992: SIGMA added
* 11-Jan-1994: Standardize
*----------------------------------------------------------------------------
      REAL*8    PI,RAD,AMP,PHASE,A11,A12,A22,B1,B2,D,C,S,XMEAN,XRMS,WMAX
      INTEGER*4 ISO,NOUT,JMAX,J
      PARAMETER (PI=3.14159265358979D0,RAD=PI/180D0)
      ISO=0
*
* Start isolating spectral lines one by one. First compute mean and variance.
*
10    CALL STATIS (N,X,XMEAN,XRMS,SIGMA(ISO))
      IF (VERBOSE .AND. ISO.GT.0) THEN
         AMP=SQRT(A(ISO)**2+B(ISO)**2)
         PHASE=ATAN2(B(ISO),A(ISO))/RAD
         WRITE (6,*) ISO,1/F(ISO),AMP,PHASE,PROB(ISO),SIGMA(ISO)
      ENDIF
*
* Check if either standard deviation criterion is reached or maximum
* number of iterations is performed.
*
      IF (SIGMA(ISO)/SIGMA(0).LT.RESID .OR. ISO.EQ.M) GOTO 300
*
* Compute periodogram
*
      CALL SPFPER(T,X,N,XMEAN,SIGMA(ISO),OFAC,HIFAC,
     .            FL,PL,L,NOUT,JMAX,PROB(ISO+1))
*
* Check the significance of the largest peak value.
*
      IF (PROB(ISO+1).GT.FALSE) GOTO 300
*
* Fit isolated frequency through data and subtract signal.
*
      ISO=ISO+1
      F(ISO)=FL(JMAX+1)
      WMAX=2*PI*F(ISO)
      IF (WMAX.EQ.0) THEN
         A(ISO)=XMEAN
         B(ISO)=0
         DO J=1,N
            X(J)=X(J)-XMEAN
         ENDDO
      ELSE
         A11=0
         A12=0
         A22=0
         B1=0
         B2=0
         DO J=1,N
            C=DCOS(WMAX*T(J))
            S=DSIN(WMAX*T(J))
            A11=A11+C*C
            A12=A12+S*C
            A22=A22+S*S
            B1=B1+C*(X(J)-XMEAN)
            B2=B2+S*(X(J)-XMEAN)
         ENDDO
         D=A11*A22-A12*A12
         A(ISO)=(+A22*B1-A12*B2)/D
         B(ISO)=(-A12*B1+A11*B2)/D
         DO J=1,N
            X(J)=X(J)-A(ISO)*DCOS(WMAX*T(J))-B(ISO)*DSIN(WMAX*T(J))
         ENDDO
      ENDIF
      GOTO 10

  300 M=ISO
      L=NOUT+1
      DO J=1,L
         PL(J)=10*LOG10(PL(J))
      ENDDO
      END
