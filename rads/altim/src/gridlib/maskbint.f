**MASKBINT -- Interogate mask in buffer
*+
      FUNCTION MASKBINT (POINTER, X, Y)
      INTEGER  MASKBINT, POINTER
      REAL*8   X, Y

* This function determines whether the mask bit at location (X,Y)
* is set or not. The mask is loaded by MASKBUFF and pointed to by
* POINTER.
*
* The location at which the mask is interogated is given by X and
* Y. These arguments are given in "world coordinates"; in other words, X
* must be between XMIN and XMAX, the coordinates of the left- and
* right-most mask point. Something similar holds for the Y-coordinate.
*
* Upon exit, the function value MASKBINT will be either 1 or 0,
* depending on the value of the mask at the pixel in which (X, Y)
* resides.
*
* When X and/or Y are out of the limits of the mask a value -1 is
* returned.
*
* Input arguments:
*  POINTER  : Pointer to the mask structure as returned by GRIDBUFF
*  X, Y     : X- and Y-coordinate of the point to be interpolated
*
* Output argument:
*  MASKBINT : Mask bit at the location (X, Y)
*             0 = bit was not set
*             1 = bit was set
*            -1 = error occurred (e.g. X or Y out of range)
*-
* 26-Jul-2006 - Avoiding use of %val
* 22-Jul-2005 - Use memory allocation
* 28-Jan-2002 - Created from GRIDBINT
*-----------------------------------------------------------------------
      integer*4 maskbints,pntr2,memloc
      include "gridbuff.inc"

* Get information about the grid

      call memget(pointer,mhead,head)
      pntr2=pointer+mhead-memloc(tmp_b)

      if (x.ge.xmin .and. x.le.xmax .and. y.ge.ymin .and. y.le.ymax) then
         
* Do the interogation

         maskbint=maskbints(tmp_i(pntr2/4),nx,ny,	! buffer, mask x,y-size,
     |   (x-xmin)/dx,(ymax-y)/dy)			! mask x,y-coord.
      else

* If the point to be interpolated is outside the mask area, or NaN, return -1

         maskbint=-1
      endif

      end

      function maskbints(z,nx,ny,x,y)
      integer maskbints
      real*8  x,y
      integer nx,ny,kx,ky,z(0:*),i
      logical btest

      kx=min(nx-1,int(x))
      ky=min(ny-1,int(y))
      i=ky*nx+kx
      if (btest(z(i/32),31-mod(i,32))) then
         maskbints=1
      else
         maskbints=0
      endif
      end
