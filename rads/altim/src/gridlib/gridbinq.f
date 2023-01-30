**GRIDBINQ -- Look-up value in buffered grid
*+
      FUNCTION GRIDBINQ (POINTER, X, Y)
      INTEGER*4	POINTER
      REAL*8	GRIDBINQ, X, Y

* This function looks up a single value in a buffered grid that
* was previously loaded using GRIDBUFF. No interpolation is
* performed, the value at the closest grid point is returned.
*
* The location at which the grid is to be queries is given by X and
* Y. These arguments are given in "world coordinates"; in other words, X
* must be between XMIN and XMAX, the coordinates of the left- and
* right-most grid point. Something similar holds for the Y-coordinate.
*
* Upon exit, the function value GRIDBINQ will be the value
* of the grid at a grid node closest to (X, Y). When that point is
* undetermined or when X and/or Y are out of limits, GRIDBINQ
* will return a NaN value.
*
* Input arguments:
*  POINTER  : Pointer to the grid structure as returned by GRIDBUFF
*  X, Y     : X- and Y-coordinate of the point to be interpolated
*
* Output argument:
*  GRIDBINQ : Value at the location (X, Y)
*-
* 26-Jul-2006 - Avoiding use of %val
*  9-Nov-2005 - Input (X, Y), not (KX, KY)
* 15-Aug-2005 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      real*8	value
      integer*4 jx,jy,pntr2,memloc
      logical	isnan
      include "gridbuff.inc"
      include "nan.inc"

* Get information about the grid

      call memget(pointer,mhead,head)
      pntr2=pointer+mhead-memloc(tmp_b)

* If the point is beyond X or Y range, return NaN

      if (x.lt.xmin .or. x.gt.xmax .or. y.lt.ymin .or. y.gt.ymax) then
         gridbinq=nan
	 return
      endif

* Determine JX (0->NX-1) and JY (0->NY-1) of the closest point

      jx=nint((x-xmin)/dx)
      jy=nint((y-ymin)/dy)

* Lookup the value

      if (ntype.eq.1) then	! BYTE
	 value=tmp_b(pntr2+jy*nx+jx)
      else if (ntype.eq.3) then	! SHORT
	 value=tmp_s(pntr2/2+jy*nx+jx)
      else if (ntype.eq.4) then	! INT
	 value=tmp_i(pntr2/4+jy*nx+jx)
      else if (ntype.eq.5) then	! FLOAT
	 value=tmp_f(pntr2/4+jy*nx+jx)
      else			! DOUBLE
	 value=tmp_d(pntr2/8+jy*nx+jx)
      endif

* Check against missing value and scale

      if (value.eq.znan .or. isnan(value)) then
         gridbinq=nan
      else
      	 gridbinq=value*dz+z0
      endif
      end
