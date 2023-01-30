**GRIDBINT -- Bi-linear interpolation of buffered grid
*+
      FUNCTION GRIDBINT (POINTER, X, Y)
      INTEGER*4	POINTER
      REAL*8	GRIDBINT, X, Y

* This function interpolates a buffered grid that was previously loaded
* using GRIDBUFF. Bi-linear interpolation is used whenever possible.
*
* The location at which the grid is to be interpolated is given by X and
* Y. These arguments are given in "world coordinates"; in other words, X
* must be between XMIN and XMAX, the coordinates of the left- and
* right-most grid point. Something similar holds for the Y-coordinate.
*
* Upon exit, the function value GRIDBINT will be the interpolated value
* of the grid at the location (X, Y). When (X, Y) points directly to a
* grid point, the value at this grid point is returned, otherwise
* bi-linear interpolation is performed between the four grid points
* surrounding (X, Y). When one of the grid points is not determined, a
* triangular interpolation is conducted.
*
* When X and/or Y are out of the limits of the grid, or when (X, Y) is
* close to an undetermined grid point, GRIDBINT will return a NaN value.
*
* Input arguments:
*  POINTER  : Pointer to the grid structure as returned by GRIDBUFF
*  X, Y     : X- and Y-coordinate of the point to be interpolated
*
* Output argument:
*  GRIDBINT : Interpolated value at the location (X, Y)
*-
* 27-Jul-2006 - Avoiding use of %val
* 18-Jan-2006 - Return NaN on NaN input
* 22-Jul-2005 - Use allocated memory
* 18-Feb-2003 - Allow for 4-byte grids
* 21-Nov-2000 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      real*8	xj,yj,value(4),weight(4),wtot,vtot
      integer*4 i,jx,jy,pntr2,memloc
      logical	isnan
      include "gridbuff.inc"
      include "nan.inc"

* Get information about the grid

      call memget(pointer,mhead,head)
      pntr2=pointer+mhead-memloc(tmp_b)

* If X or Y are NaN or beyond allowed range, return NaN

      if (x.ge.xmin.and.x.le.xmax .and. y.ge.ymin.and.y.le.ymax) then
      else
         gridbint=nan
	 return
      endif

* Determine JX (0->NX-1) and JY (0->NY-1) of lower left corner

      xj=(x-xmin)/dx
      yj=(y-ymin)/dy
      jx=min(nx-2,int(xj))
      jy=min(ny-2,int(yj))
      xj=xj-jx
      yj=yj-jy

* Lookup 4 corners

      if (ntype.eq.1) then	! BYTE
	 pntr2=pntr2+jy*nx+jx
	 value(1)=tmp_b(pntr2)
	 value(2)=tmp_b(pntr2+1)
	 value(3)=tmp_b(pntr2+nx)
	 value(4)=tmp_b(pntr2+nx+1)
      else if (ntype.eq.3) then	! SHORT
	 pntr2=pntr2/2+jy*nx+jx
	 value(1)=tmp_s(pntr2)
	 value(2)=tmp_s(pntr2+1)
	 value(3)=tmp_s(pntr2+nx)
	 value(4)=tmp_s(pntr2+nx+1)
      else if (ntype.eq.4) then	! INT
	 pntr2=pntr2/4+jy*nx+jx
	 value(1)=tmp_i(pntr2)
	 value(2)=tmp_i(pntr2+1)
	 value(3)=tmp_i(pntr2+nx)
	 value(4)=tmp_i(pntr2+nx+1)
      else if (ntype.eq.5) then	! FLOAT
	 pntr2=pntr2/4+jy*nx+jx
	 value(1)=tmp_f(pntr2)
	 value(2)=tmp_f(pntr2+1)
	 value(3)=tmp_f(pntr2+nx)
	 value(4)=tmp_f(pntr2+nx+1)
      else			! DOUBLE
	 pntr2=pntr2/8+jy*nx+jx
	 value(1)=tmp_d(pntr2)
	 value(2)=tmp_d(pntr2+1)
	 value(3)=tmp_d(pntr2+nx)
	 value(4)=tmp_d(pntr2+nx+1)
      endif

* Set corner weights

      weight(1)=(1-xj)*(1-yj)
      weight(2)=   xj *(1-yj)
      weight(3)=(1-xj)*   yj
      weight(4)=   xj *   yj

* Add up weights

      wtot=0d0
      vtot=0d0
      do i=1,4
         if (value(i).eq.znan .or. isnan(value(i))) then
	 else
	    wtot=wtot+weight(i)
	    vtot=vtot+weight(i)*value(i)
	 endif
      enddo
      gridbint=vtot/wtot*dz+z0
      end
