**F1F2 -- Convert sigma0 to wind and vv using Gourrion 2-parameter model
*+
      FUNCTION F1F2 (I, VAL, SWH)
      REAL*8    F1F2, VAL, SWH
      INTEGER*4 I
*
* This function computes wind speed (U10, referenced to 10 m above
* sea level) from altimeter backscatter coefficient (SIGMA0) and
* significant wave height (SWH). The reverse is also possible:
* converting U10 and SWH to SIGMA0.
* The parameter I selects whether the forward (I=1) or backward
* (I=2) conversion is requested.
*
* The unit for SIGMA0 is dB, for SWH is m and for U10 is m/s.
*
* The respective models F1(SIGMA0,SWH) and F2(U10,SWH) are described in
* Gourrion et al. (2002).
*
* Input arguments:
*  I   : Select direction of conversion
*        1 = Convert SIGMA0 and SWH to U10 (F1)
*        2 = Convert U10 and SWH to SIGMA0 (F2)
*  VAL : Either SIGMA0 (I=1) or U10 (I=2)
*  SWH : Significant wave height
*
* Returned value:
*  F1F2: Either U10 (I=1) or SIGMA0 (I=2)
*
* Reference:
*
* Gourrion, J., D. Vandemark, S. Bailey, B. Chapron, G. P. Gommenginger,
* P. G. Challenor, and M. A. Srokosz,
* A two-parameter wind speed algorithm for Ku-band altimeters,
* J. Atmos. Ocean. Technol., 19(12), 2030-2048, 2002.
*-
* $Log: f1f2.f,v $
* Revision 1.3  2006/07/28 22:03:53  rads
* - Removed obsolete Fortran constructions
*
* Revision 1.2  2006/01/28 19:52:08  rads
* - Return NaN on out of range input
*
* Revision 1.1  2004/08/28 20:25:41  remko
* - F1F2 added to RSSUBS suite
*
*-----------------------------------------------------------------------
      real*8	a(3)/-0.34336d0,0.10000d0,0.08725d0/
      real*8	b(3)/ 0.06909d0,0.02844d0,0.06374d0/
      ! Note: second and third element of a and b intensionally switched
      real*8	wx(2,2,2)/-33.95062d0,-11.03394d0,-3.93428d0,-0.05834d0,
     |			  -43.39541d0, -6.92550d0, 2.78612d0, 1.22293d0/
      real*8	wy(2,2)	 /  0.54012d0, 10.40481d0, 1.18281d0,-3.30096d0/
      real*8	bx(2,2)	 / 18.06378d0, -0.37228d0, 7.83459d0,-1.46489d0/
      real*8	by(2)	 / -2.28387d0,  1.13906d0/
      save	a,b,wx,wy,bx,by
      real*8	x(2),y,p(2)
      integer	j
      logical	isnan
      include "nan.inc"

      if (isnan(val) .or. isnan(swh) .or. val.gt.1d20 .or. swh.gt.1d20
     |		.or. i.lt.1 .or. i.gt.2) then
         ! If input values are NaN or VOID, or incomprehensible input, return NaN
	 f1f2=nan
	 return
      endif

* Normalise input values

      p(1)=a(i)+b(i)*val
      p(2)=a(3)+b(3)*swh

* Determine X vector

      do j=1,2
	 x(j)=wx(1,j,i)*p(1)+wx(2,j,i)*p(2)+bx(j,i)
      enddo
      do j=1,2
         x(j)=1d0/(1d0+exp(-x(j)))
      enddo

* Determine normalized output Y

      y=wy(1,i)*x(1)+wy(2,i)*x(2)+by(i)
      y=1d0/(1d0+exp(-y))

* Return unnormalized value

      f1f2=(y-a(3-i))/b(3-i)
      end
