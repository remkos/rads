**SETAREA -- Set the plot area from grid info (if needed)
*
      subroutine setarea (i,nx,ny,xg0,xg1,yg0,yg1,cell)
      integer i,nx,ny
      real xg0,xg1,yg0,yg1
      logical cell

*     call pmqwin(xw0,xw1,yw0,yw1)
      include "pim.inc"

      if (i.eq.1 .and. cell) then

* Set plotarea according to gridarea (for cells !)

         if (xw0.eq.xw1 .and. yw0.eq.yw1) then
	    xw0=xg0-(xg1-xg0)/(nx-1)/2
	    xw1=xg1+(xg1-xg0)/(nx-1)/2
	    yw0=yg0-(yg1-yg0)/(ny-1)/2
	    yw1=yg1+(yg1-yg0)/(ny-1)/2
	 endif

      else if (i.eq.1) then

* Set plotarea according to gridarea (for grids !)

         if (xw0.eq.xw1 .and. yw0.eq.yw1) then
	    xw0=xg0
	    xw1=xg1
	    yw0=yg0
	    yw1=yg1
	 endif

      else if (cell) then

* Set gridarea according to window area (for cells !)

         if (xg0.eq.xg1 .and. yg0.eq.yg1) then
	    xg0=xw0+(xw1-xw0)/(nx-1)/2
	    xg1=xw1-(xw1-xw0)/(nx-1)/2
	    yg0=yw0+(yw1-yw0)/(ny-1)/2
	    yg1=yw1-(yw1-yw0)/(ny-1)/2
	 endif

      else

* Set gridarea according to window area (for grids !)

         if (xg0.eq.xg1 .and. yg0.eq.yg1) then
	    xg0=xw0
	    xg1=xw1
	    yg0=yw0
	    yg1=yw1
	 endif

      endif
      end
