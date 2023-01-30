**HELMERT1 -- Determine Helmert transformation coefficients
*+
      FUNCTION HELMERT1 (N, XYZFROM, XYZTO, SIGMA, COEFF, WRMS)
      INTEGER HELMERT1, N
      REAL*8 XYZFROM(3,N), XYZTO(3,N), SIGMA(3,N), COEFF(7), WRMS(2)

* This routine determines the Helmert transformation coefficients by
* matching the Cartesian coordinates of N points stored in the array
* XYZFROM with those stored in XYZTO. The solution of the seven
* transformation coefficients are obtained by least squares fitting. The
* different points may be assigned different weights by means of the
* standard deviations stored in the array SIGMA. Points can be ignored
* by giving a negative SIGMA. When SIGMA(1,1)=0, all points will get equal
* weights.
*
* The seven Helmert transformation coefficients are stored in the array
* COEFF. They are:
*
* COEFF(1..3) : Translation vector (X,Y,Z) in meters
* COEFF(4)    : Scale difference
* COEFF(5..7) : Rotations around the X, Y and Z axes in radians
*
* Arguments:
*  N       (input) : Number of points to transform
*  XYZFROM (input) : XYZ-Coordinates of the points to transform
*  XYZTO   (input) : Reference XYZ-Coordinates
*  SIGMA   (input) : Standard deviation in XYZ of each of the points
*  COEFF  (output) : Helmert transformation coefficients (see above)
*  WRMS   (output) : Weighted RMS residual. Two values:
*                    WRMS(1) = a apriori, WRMS(2) = a posteriori
*  HELMERT1 (output) : Returned error value:
*                    0 = No errors, 1 = Insufficient number of points
*                    2 = Too few weighted points,
*                    3 = Error in normal matrix factorisation
*                    4 = Error in normal matrix solve
*
* Note that the position coordinates in XYZFROM and XYZTO and the
* standard deviations SIGMA all have to be in the same units as
* COEFF(1..3).
*-
*  4-Apr-2001 - Created by Remko Scharroo
*-----------------------------------------------------------------------
* The 21 partial derivatives
      real*8 dxdc(7,3)/1d0,7*0d0,1d0,7*0d0,1d0,4*0d0/
      integer npar,i,j,k,l,m,nwght
      parameter (npar=7)
      real*8 atwa(npar*(npar+1)/2),atwr(npar),rtwr
      real*8 w,r

* Initialise normal matrix (AtWA), righthand vector (AtWR) and
* weighted sum square residual (RtWR)

      do i=1,npar*(npar+1)/2
         atwa(i)=0d0
      enddo
      do i=1,npar
         atwr(i)=0d0
      enddo
      rtwr=0d0
      nwght=0

* Check if sufficient number of points are given

      if (n.lt.3) then
         write (*,*) 'HELMERT1: insufficient number of points'
	 helmert1=1
	 return
      endif

* Loop through all the points

      do i=1,n

* Determine non-zero and non-unity partials

	 dxdc(4,1)= xyzto(1,i)
	 dxdc(6,1)= xyzto(3,i)
	 dxdc(7,1)=-xyzto(2,i)
	 dxdc(4,2)= xyzto(2,i)
	 dxdc(5,2)=-xyzto(3,i)
	 dxdc(7,2)= xyzto(1,i)
	 dxdc(4,3)= xyzto(3,i)
	 dxdc(5,3)= xyzto(2,i)
	 dxdc(6,3)=-xyzto(1,i)

* Update AtWA, AtWR, and RtWR

         do j=1,3	! Cycle through 3 coordinates
	    if (sigma(1,1).eq.0d0) then
	       w=1	! When sigma(1,1)=0 assign all point same weight
            else if (sigma(j,i).lt.0d0) then
	       goto 200	! Ignore coordinates with sigma(j,i)<0
	    else
	       w=1/sigma(j,i)**2
	    endif
	    k=0
	    r=xyzto(j,i)-xyzfrom(j,i)	! Residual in jth coordinate
	    do l=1,npar
	       do m=1,l
	          k=k+1
	          atwa(k)=atwa(k)+dxdc(l,j)*w*dxdc(m,j)
	       enddo
	       atwr(l)=atwr(l)+dxdc(l,j)*w*r
	    enddo
	    rtwr=rtwr+r*w*r
	    nwght=nwght+1	! Count weighted coordinates
200	    continue
        enddo
      enddo

* Check if there are enough weighted measurements

      if (nwght.lt.7) then
	 write (*,*) 'HELMERT1: too few weighted points'
	 helmert1=2
	 return
      endif

* Solve least-squares problem
* Step 1: Factorise the normal matrix (AtWA)

      call dpptrf('U',npar,atwa,i)
      if (i.ne.0) then
         write (*,*) 'HELMERT1: error in matrix factorisation'
	 helmert1=3
	 return
      endif

* Step 2: Store the right hand vector (AtWR) for later use and solve
*         x from AtWA x = AtWR

      do i=1,npar
         coeff(i)=atwr(i)
      enddo
      call dpptrs('U',npar,1,atwa,coeff,npar,i)
      if (i.ne.0) then
	 write (*,*) 'HELMERT1: error in matrix solve'
	 helmert1=4
	 return
      endif

* Compute a priori and posteriori weighted RMS

      wrms(1)=sqrt(rtwr/nwght)
      do i=1,npar
         rtwr=rtwr-coeff(i)*atwr(i)
      enddo
      wrms(2)=sqrt(rtwr/nwght)
      helmert1=0
      end
