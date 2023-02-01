**RESFITEDIT -- Edit out points beyond a certain distance from the fit
*+
      SUBROUTINE RESFITEDIT (N, X, Y, SIGMA, A, B, FIT, EDIT)
      INTEGER N
      REAL*8 X(N), Y(N), SIGMA(N), A, B, FIT, EDIT

* This routine edits out points Y (function of X) that are beyond a
* certain threshold away from a predetermined fit of the points in
* an X-Y diagram.
*
* The edit criterion is |(A+B*X)-Y| > EDIT*FIT where FIT is the RMS of
* fit of all (selected points, i.e. SIGMA>0) and EDIT is an edit multiplier.
*
* Points that are edited out will be given a negative SIGMA. But ALL points
* are tested against the threshold, also those that have been edited
* before. The magnitude of SIGMA itself is not used.
*
* Arguments:
*  N     : Number of data points
*  X     : Array of X values of the data points
*  Y     : Array of Y values of the data points
*  SIGMA : Array of Y error variance of the data points
*  A     : Bias of the fit
*  B     : Tilt of the fit
*  FIT   : RMS of fit
*-
* 25-Aug-1999 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer i

      do i=1,N
         if (abs(a+b*x(i)-y(i)).gt.edit*fit) then
	    sigma(i)=-abs(sigma(i))
	 else
	    sigma(i)=+abs(sigma(i))
	 endif
      enddo
      end
