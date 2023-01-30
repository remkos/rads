**PIXINTER -- Interpolate grid to pixel resolution
*
      subroutine pixinter(a,nx,ny,x0,x1,y0,y1,b,px,py)
      integer nx,ny,px,py
      real*4 a(nx,ny),b(px,py)
      integer ix,iy
      real gr1int4,x,y,xfact1,xfact2,yfact1,yfact2
      real xblc,xtrc,yblc,ytrc,x0,x1,y0,y1

* Note: pixels are cells !
* xfact1,yfact1: conversion from pixels to window coordinates
* xfact2,yfact2: conversion from geodetic coordinates to grid coordinates

      call pgqwin(xblc,xtrc,yblc,ytrc)
      xfact1=(xtrc-xblc)/px
      yfact1=(ytrc-yblc)/py
      xfact2=(nx-1)/(x1-x0)
      yfact2=(ny-1)/(y1-y0)
      do iy=1,py
         do ix=1,px
	    x=xblc+(ix-0.5)*xfact1
	    y=yblc+(iy-0.5)*yfact1
            call pmcinv(1,x,y)
            x=(x-x0)*xfact2+1
            y=(y-y0)*yfact2+1
            b(ix,iy)=gr1int4(a,x,y,nx,ny,nx)
	 enddo
      enddo
      end
