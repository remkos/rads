**RESFIT -- Make a fit through a set of X and Y values
*+
      SUBROUTINE RESFIT (NOBS, X, Y, SIGMA, IOBS, YMEAN, YRMS, ERMS,
     |		A, B, FIT)
      INTEGER NOBS
      REAL*8 X(NOBS), Y(NOBS), SIGMA(NOBS)
      INTEGER IOBS
      REAL*8 YMEAN, YRMS, ERMS, A, B, FIT
*-
* 24-Aug-1999 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      real*8 xmean, xrms, sumxy, uxx, uxy, uyy
      integer i

      IOBS=0
      xmean=0d0
      YMEAN=0d0
      sumxy=0d0
      xrms=0d0
      YRMS=0d0
      ERMS=0d0

* Compute the statistics of the observation residuals for this pass.
* Computed are: mean, rms, fit (linear regression in range-rangerate
* diagram), and rms of fit.
* Only valid measurements are used.

      do i=1,nobs
	 if (sigma(i).gt.0d0) then
	    IOBS=IOBS+1
            xmean=xmean+X(i)
            YMEAN=YMEAN+Y(i)
            sumxy=sumxy+X(i)*Y(i)
            xrms =xrms +X(i)*X(i)
            YRMS =YRMS +Y(i)*Y(i)
	 else
	    ERMS =ERMS +Y(i)*Y(i)
	 endif
      enddo

* Statistics are determined, now compute parameters of the fit:
* a (bias), b (tilt), fit (rms of fit)

      if (IOBS.le.1) then
	 A=YMEAN
	 B=0
	 fit=0
      else
         uxx=IOBS*xrms -xmean*xmean
         uxy=IOBS*sumxy-xmean*YMEAN
         uyy=IOBS*YRMS -YMEAN*YMEAN
         B=uxy/uxx
         A=(YMEAN-B*xmean)/IOBS
   
	 FIT=0
	 do i=1,nobs
	    if (sigma(i).gt.0d0) then
	       FIT=FIT+(A+B*X(i)-Y(i))**2
	    endif
	 enddo
      endif
      if (IOBS.ge.1) then
	 FIT=sqrt(FIT/IOBS)
         YMEAN=YMEAN/IOBS
         YRMS=sqrt(YRMS/IOBS)
      endif
      if (IOBS.NE.NOBS) ERMS=sqrt(ERMS/(NOBS-IOBS))

      end
