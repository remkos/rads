**REGRES -- Regression analysis
*+
      SUBROUTINE REGRES (N, X, Y, A, B, R)
      INTEGER*4 N
      REAL*8  X(N), Y(N), A, B, R
*
* This subroutine computes the best fitting straight line trough a number
* of points with coordinates (X,Y). Upon return, A and B will be the
* coefficients of the line
*
*     Y = A + B * X
*
* that is the best fit through the data points. R is the correlation
* between the fit and the data.
*
* Arguments:
*  N  (input) : Number of data points
*  X  (input) : Array of X coordinates of the data points
*  Y  (input) : Array of Y coordinates of the data points
*  A (output) : Coefficient of linear fit (Y = A + B * X)
*  B (output) : Coefficient of linear fit (Y = A + B * X)
*  R (output) : Correlation
*-
*  2-Jul-1993 - New manual
* 13-Apr-2005 - Skip NaNs
*-----------------------------------------------------------------------
      integer*4 i,m
      real*8	sumx,sumy,sumxx,sumxy,sumyy,uxx,uxy,uyy
      logical	isnan
      sumx=0d0
      sumy=0d0
      sumxy=0d0
      sumxx=0d0
      sumyy=0d0
      m=0
      do i=1,n
	 if (.not.isnan(x(i)).and..not.isnan(y(i))) then
	    m=m+1
            sumx =sumx +x(i)
            sumy =sumy +y(i)
            sumxy=sumxy+x(i)*y(i)
            sumxx=sumxx+x(i)*x(i)
            sumyy=sumyy+y(i)*y(i)
         endif
      enddo
      uxx=m*sumxx-sumx*sumx
      uxy=m*sumxy-sumx*sumy
      uyy=m*sumyy-sumy*sumy
      b=uxy/uxx
      a=(sumy-b*sumx)/m
      r=uxy/sqrt(uxx*uyy)
      end
