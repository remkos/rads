**STATIS -- Compute statistics of a series
*+
      SUBROUTINE STATIS (N, X, XMEAN, XRMS, XSIGMA)
      INTEGER*4 N
      REAL*8  X(N), XMEAN, XRMS, XSIGMA
*
* This routine computes the average (XMEAN), root-mean-squate (XRMS) and
* standard deviation (XSIGMA) of a series of N values (X).
*
* Arguments:
*  N      (input) : Number of values in the series.
*  X      (input) : Series of values.
*  XMEAN (output) : Average of the series.
*  XRMS  (output) : RMS of the series.
*  XSIGMA (output) : Standard deviation of the series.
*-
* 18-Dec-1992 - Created: Remko Scharroo, DUT/SOM (c)
* 11-Jan-1994 - Standardize.
*-----------------------------------------------------------------------
      REAL*8 SUM, SUM2
      INTEGER*4 I

      SUM=0
      SUM2=0

      DO I=1,N
	 SUM=SUM+X(I)
	 SUM2=SUM2+X(I)**2
      ENDDO
      XMEAN=SUM/N
      XRMS=SQRT(SUM2/N)
      XSIGMA=SQRT((SUM2-SUM**2/N)/(N-1))
      END
