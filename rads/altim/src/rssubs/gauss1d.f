**GAUSS1D -- One-dimensional Gaussian filter
*+
      SUBROUTINE GAUSS1D (N, X, SIGX, HORX, Z, SIGZ, ZF)
      INTEGER*4 N
      REAL*8    X(N), SIGX, HORX, Z(N), SIGZ(N), ZF(N)

* This subroutine filters a set of N data points aligned along a (semi-)
* one-dimensional coordinate. The independent variable is X, the data
* values are N. All data points should be ordered in ascending or
* descending order of X.
*
* The data are weighted by their respective variances SIGZ,
* and the spacial scale SIGX.
* - Data with negative SIGZ are not used, but the filtered value is return.
* - When SIGZ(1)=0, all data points are asigned an equal variance: SIGZ(i)=1
*
* The horizon HORX determines the window to
* which the filtering is applied (usually HORX = 2.5 * SIGX).
* The (low-pass) filtered data are returned as ZF.
*
* The Guassian filtering is formulated as follows:
*
* ZF(i) = sum [ w(i,j) * Z(i) ] / sum w(i,j)
*          j                       j
*
* where w(i,j) is the weight contribution of point j to point i, given by
*
* w(i,j) = exp [ -(X(i)-X(j))**2 / SIGX**2 ] / SIGZ(j)
*
* ZF(i) is set to 1d30 when the value can not be determined.
*
* Arguments:
* N    (input) : Number of data points
* X    (input) : Independent variable
* SIGX (input) : Distance weighting sigma
* HORX (input) : Horizon (best value: 2.5*SIGX)
* Z    (input) : Data values
* SIGZ (input) : Data variances. When SIGZ(i)<0, point i is not used in
*	         the filtering, but ZF(i) will be defined.
*                When SIGZ(1)=0, SIGZ(i)=1 will be used for all i.
* ZF  (output) : Low-pass filtered data
*-
* 26-Oct-95 : Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4	i,j,j0
      real*8	dx2,horx2,sigx2,wgt,sum,w

      j0=1
      horx2=horx**2
      sigx2=sigx**2

      do i=1,n
	 sum=0d0
	 wgt=0d0
	 do j=j0,n
	    dx2=(x(i)-x(j))**2
	    if (dx2.gt.horx2) then
	       if (j.lt.i) then
	          j0=j+1
	       else
	          goto 100
	       endif
	    else if (sigz(1).eq.0) then
	       w=exp(-dx2/sigx2)
	       wgt=wgt+w
	       sum=sum+w*z(j)
	    else if (sigz(j).gt.0) then
	       w=exp(-dx2/sigx2)/sigz(j)
	       wgt=wgt+w
	       sum=sum+w*z(j)
	    endif
	 enddo
100	 continue
	 if (wgt.eq.0d0) then
	    zf(i)=1d30
	 else
	    zf(i)=sum/wgt
	 endif
      enddo
      end
