**GRDINTER -- Interpolate grid to lower resolution
*
      subroutine grdinter(a,nx,ny,b,px,py)
      integer nx,ny,px,py
      real*4 a(nx,ny),b(px,py)
      integer ix,iy
      real gr1int4,x,y,xfact,yfact

      xfact=real(nx-1)/real(px-1)
      yfact=real(ny-1)/real(py-1)
      do iy=1,py
	 y=(iy-1)*yfact+1
         do ix=1,px
	    x=(ix-1)*xfact+1
            b(ix,iy)=gr1int4(a,x,y,nx,ny,nx)
	 enddo
      enddo
      end
