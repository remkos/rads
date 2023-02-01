**DYNTOPO -- Compute velocity from dynamic topography
*+
      SUBROUTINE DYNTOPO (Z, VX, VY, NX, NY, MX, X0, X1, Y0, Y1)
      INTEGER NX, NY, MX
      REAL Z(MX,*), VX(MX,*), VY(MX,*), X0, X1, Y0, Y1
*
* Determine velocity field components VX (due West) and VY (due North)
* based on the dynamic topography height Z. The velocity grids and
* the dynamic topography grid have the same resolution and are for the
* same area.
*
* Arguments:
*  Z       (input): grid of dynamic height (unit: metres)
*  VX, VY (output): grids of X and Y component of the velocity field
*                   (unit: metres per second)
*  NX, NY  (input): dimension of grids
*  MX      (input): first dimension of arrays Z, VX, and VY
*  X0, X1  (input): longitude boundaries (deg)
*  Y0, Y1  (input): latitude boundaries (deg)
*-
*  3-Jan-1995 - Created by Remko Scharroo
*  5-Jan-1996 - New manual.
*-----------------------------------------------------------------------
      integer kx,ky
      real y,dx,dy,sx,sy
      real rad,dynmul

      real rotate,am,gm,dynfac
      parameter (rotate=7.292115855e-5,am=6371000e0)
      parameter (gm=9.8062,dynfac=gm/(2*rotate))
      
      rad=4*atan(1e0)/180

      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)

      call griddy(z,vx,nx,ny,nx)
      call griddx(z,vy,nx,ny,nx)

      sy=am*dy*rad
      do ky=1,ny
         y=y0+(ky-1)*dy
         sx=am*dx*rad*cos(y*rad)
         dynmul=dynfac/sin(y*rad)
         do kx=1,nx
	    vx(kx,ky)=-dynmul*vx(kx,ky)/sy
	    vy(kx,ky)= dynmul*vy(kx,ky)/sx
         enddo
      enddo

      end
