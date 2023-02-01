**VELOCITY -- Plot vectors according to velocity grids
*+
      SUBROUTINE VELOCITY (VX, VY, NX, NY, MX, MY, X0, X1, Y0, Y1,
     &			   V0, V1, FACTOR)
      INTEGER NX, NY, MX, MY
      REAL VX(NX,*), VY(NX,*), X0, X1, Y0, Y1, V0, V1, FACTOR
*
* Plot velocity field components VX and VY as arrows in the
* plot area. The arrows are sized depending on the velocity and aligned
* depending on the X and Y component of the velocity fields.
* The size of the arrow head is determined by the current character height.
*
* Arguments (input):
*  VX, VY : grids of X and Y component of the velocity field
*  NX, NY : dimension of the grids
*  MX, MY : step size for plotting vectors
*  X0, X1 : longitude boundaries of the grids
*  Y0, Y1 : latitude boundaries of the grids
*  V0, V1 : cutoff values of velocity for plotting vectors
*  FACTOR : factor transforming velocity (m/s) to vector size (m),
*           e.g. FACTOR=1e-2 means 1 m/s is plotted as 1 cm.
*
      integer kx,ky
      real x,y,dx,dy,angle,fx,fy,v,c,s,length
      real rad
      real xw0,xw1,yw0,yw1
      real xp0,xp1,yp0,yp1
      real cutoff
      parameter (cutoff=-10.0)

      call pgqwin(xw0,xw1,yw0,yw1)
      call pgqvp(2,xp0,xp1,yp0,yp1)
      
      rad=4*atan(1e0)/180
      
      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)
      fx=(xw1-xw0)/(xp1-xp0)
      fy=(yw1-yw0)/(yp1-yp0)

      call pgsah(2,60.0,1.0)

      call pgbbuf
      do ky=1,ny,my
         y=y0+(ky-1)*dy
	 if (abs(y).ge.cutoff) then
         do kx=1,nx,mx
	    v=sqrt(vx(kx,ky)**2+vy(kx,ky)**2)
	    if (v.ge.v0 .and. v.le.v1) then
               x=x0+(kx-1)*dx
               y=y0+(ky-1)*dy
               angle=atan2(vx(kx,ky),vy(kx,ky))/rad
	       call pmcvec(1,x,y,angle)
               c=cos(angle*rad)
               s=sin(angle*rad)
	       length=v*factor*1e3/2
	       call pgarro(
     |		  x-c*length*fx,
     |		  y-s*length*fy,
     |		  x+c*length*fx,
     |		  y+s*length*fy)
	    endif
         enddo
	 endif
      enddo
      call pgebuf

      end
