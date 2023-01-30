**SPFPER -- Compute Fast Lomb-Scargle normalized periodogram.
*+
      SUBROUTINE SPFPER (X, Y, N, YMEAN, YSIGMA, OFAC, HIFAC,
     |                   WK1, WK2, NWK, NOUT, JMAX, PROB)
      INTEGER*4 N, NWK, NOUT, JMAX
      REAL*8    X(N), Y(N), YMEAN, YSIGMA, OFAC, HIFAC,
     |          WK1(0:NWK-1), WK2(0:NWK-1), PROB
*
* Given N data points with abscissas X (which need not be equally spaced,
* but must be sequential) and ordinates Y (with mean YMEAN and standard dev
* YSIGMA), and given a desired oversampling factor OFAC (a typical value
* being 4 or larger), this routine computes the normalized periodogram
* of the data points using the Lomb-Scargle method (see reference below).
*
* The routine fills array WK1 with a sequence of NOUT increasing frequencies
* (not angular frequencies but in units of 1 over the unit of X) up to HIFAC
* times the `average' Nyquist frequency, and fills array WK2 with the values
* of the Lomb-Scargle normalized periodogram at those frequencies. The
* arrays X and Y are not altered.
*
* NWK, the dimension of WK1 and WK2, must be large enough for intermediate
* work space, or an error (pause) results. The routine also returns JMAX,
* such that WK2(JMAX) is the maximum element in WK2, and PROB, an estimate
* of the significance of that maximum against the hypothesis of random noise.
* A small value of PROB indicates that a significant periodic signal is
* present.
*
* In order to infer the amplitudes AMP from the normalized periodogram WK2,
* compute AMP(K) = 2 * YSIGMA * SQRT(WK2(K)/N).
*
* Arguments:
*   X      (input): abscissas of the data points (e.g. time).
*   Y      (input): N real data points, need not be equally spaced.
*   N      (input): Number of data points.
*   YMEAN  (input): Average of the N data values.
*   YSIGMA (input): Standard deviation of the N data values.
*   OFAC   (input): Oversampling factor.
*   HIFAC  (input): Determines the highest frequency in the periodogram (in
*                   average Nyquist units).
*   WK1   (output): Frequencies (NOT angular, in 1/units of X)
*   WK2   (output): Normalized periodogram.
*   NWK    (input): Size of the working spaces as defined in the calling
*                   (sub)program.
*   NOUT  (output): Number of peaks in periodogram.
*   JMAX  (output): Index of highest peak in the periodogram.
*   PROB  (output): Significance of the highest peak in the periodogram.
*
* Ref: Press W.H. and George B. Rybicki.  Fast algorithm for spectral analysis
* of unevenly sampled data, The Astronomical Journal, 338, 227-280, 1989.
*-
* 11-Mar-1991: Created - Remko Scharroo
* 10-Dec-1992: Bug fixed (XMIN -> X(1)) [RS]
* 11-Mar-1993: New manual.
* 11-Jan-1994: Standardize.
*-----------------------------------------------------------------------
      INTEGER MACC,J,JPNT,K,NPNT
      REAL*8  EFFM,CTERM,STERM,WTAU,CWTAU,SWTAU,C2W,C2WTJ,S2WTJ,HCWTJ,
     |        HSWTJ,PMAX,XDIF,FAC,FNPNT,CK,CKK,DF,EXPY
      PARAMETER (MACC=4)
      NOUT=.5D0*OFAC*HIFAC*N
      JPNT=NINT(OFAC*HIFAC*N*MACC)
      NPNT=64
   10 IF (NPNT.LT.JPNT) THEN
         NPNT=NPNT*2
         GOTO 10
      ENDIF
      IF (NPNT+1.GT.NWK-1) THEN
         WRITE (0,600) NPNT+2
         STOP
      ENDIF
*
* Set workspaces to zero.
*
      DO J=0,NPNT+1
         WK1(J)=0
         WK2(J)=0
      ENDDO
*
* `Extirpolate' the data into the workspaces and take their reverse FFTs.
*
      XDIF=X(N)-X(1)
      FAC=NPNT/XDIF/OFAC
      FNPNT=NPNT
      DO J=1,N
         CK=1+DMOD((X(J)-X(1))*FAC,FNPNT)
         CKK=1+DMOD(2*CK-2,FNPNT)
         CALL SPREAD1(Y(J)-YMEAN,WK1,NPNT,CK,MACC)
         CALL SPREAD1(1D0,WK2,NPNT,CKK,MACC)
      ENDDO
      CALL SPFFTR(WK1,NPNT,-1)
      CALL SPFFTR(WK2,NPNT,-1)
*
* Compute the Lomb-Scargle value for each frequency.
*
      DF=(N-1)/XDIF/OFAC/N
      WK1(0)=0
      WK2(0)=(YMEAN/YSIGMA)**2
      JMAX=0
      PMAX=WK2(0)
      DO J=1,NOUT
         K=2*J
         HCWTJ=WK1(K)
         HSWTJ=-WK1(K+1)
         C2WTJ=WK2(K)
         S2WTJ=-WK2(K+1)
         WTAU=DATAN2(S2WTJ,C2WTJ)
         C2W=C2WTJ*DCOS(2*WTAU)+S2WTJ*DSIN(2*WTAU)
         CWTAU=DCOS(WTAU)
         SWTAU=DSIN(WTAU)
         CTERM=(HCWTJ*CWTAU+HSWTJ*SWTAU)**2/(N+C2W)
         STERM=(HSWTJ*CWTAU-HCWTJ*SWTAU)**2/(N-C2W)
         WK1(J)=J*DF
         WK2(J)=(CTERM+STERM)/YSIGMA**2
         IF (WK2(J).GT.PMAX) THEN
            PMAX=WK2(J)
            JMAX=J
         ENDIF
      ENDDO
*
* Estimate the significance of the largest peak value.
*
      EXPY=DEXP(-PMAX)
      EFFM=2*NOUT/OFAC
      PROB=EFFM*EXPY
      IF (PROB.GT.1D-2) PROB=1-(1-EXPY)**EFFM
*
* Formats
*
  600 FORMAT ('spfper: Working space too small, need at least',i9)
      END


      SUBROUTINE SPREAD1 (Y, YY, N, X, M)
      INTEGER IHI,ILO,IX,J,M,N,NDEN
      REAL*8  X,Y,FAC
      REAL*8  YY(N)
*
* Given an array YY of length N, extirpolate (spread) a value Y into M actual
* array elements that best approximate the `fictional' (i.e. possibly non-
* integer) array element number X. The weights used are coefficients of the
* Lagrange interpolation polynomial.
*
      INTEGER NFAC(10) /1,1,2,6,24,120,720,5040,40320,362880/
      IF (M.GT.10) STOP 'spread: Factorial table too small'
      IX=X
      IF (X.EQ.DFLOAT(IX)) THEN
         YY(IX)=YY(IX)+Y
      ELSE
         ILO=MIN(MAX(INT(X-.5D0*M+1),1),N-M+1)
         IHI=ILO+M-1
         NDEN=NFAC(M)
         FAC=X-ILO
         DO J=ILO+1,IHI
            FAC=FAC*(X-J)
         ENDDO
         YY(IHI)=YY(IHI)+Y*FAC/NDEN/(X-IHI)
         DO J=IHI-1,ILO,-1
            NDEN=NDEN/(J+1-ILO)*(J-IHI)
            YY(J)=YY(J)+Y*FAC/NDEN/(X-J)
         ENDDO
      ENDIF
      END
