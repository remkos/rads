**GRIDBSPL -- Bi-cubic spline interpolation of buffered grid
*+
      FUNCTION GRIDBSPL (POINTER, X, Y)
      INTEGER*4	POINTER
      REAL*8	GRIDBSPL, X, Y

* This function interpolates a buffered grid that was previously loaded
* using GRIDBUFF. Bi-cubic spline interpolation is used whenever possible.
*
* The location at which the grid is to be interpolated is given by X and
* Y. These arguments are given in "world coordinates"; in other words, X
* must be between XMIN and XMAX, the coordinates of the left- and
* right-most grid point. Something similar holds for the Y-coordinate.
*
* Upon exit, the function value GRIDBSPL will be the interpolated value
* of the grid at the location (X, Y). When (X, Y) points directly to a
* grid point, the value at this grid point is returned, otherwise
* bi-cubic spline interpolation is performed between the 6-by-6 grid
* points surrounding (X, Y).
*
* This routine DOES NOT take into account invalid values. Hence, it
* assumes that all values in the 6-by-6 subgrid are valid. This will
* work for most geoid and mean sea surface grids.
*
* When X and/or Y are out of the limits of the grid, or too close to the
* boundary in order to perform bi-cubic spline interpolation,
* GRIDBSPL will return a NaN value.
*
* Input arguments:
*  POINTER : Pointer to the grid structure as returned by GRIDBUFF
*  X, Y    : X- and Y-coordinate of the point to be interpolated
*
* Output argument:
*  GRIDBSPL : Interpolated value at the location (X, Y)
*-
* 27-Jul-2006 - Avoiding use of %val
* 22-Jul-2005 - Use allocated memory
* 18-Feb-2003 - Created from GRIDBINT and Numerical Recipes
*-----------------------------------------------------------------------
      integer*4	jx,jy,nwin,nh,kx,ky,pntr2,memloc
      parameter (nwin=6,nh=nwin/2)
      real*8	xj,yj,z(nwin,nwin),zy(nwin),w(nwin),u(nwin),gridbsplu
      include "gridbuff.inc"
      include "nan.inc"

* Get information about the grid

      call memget(pointer,mhead,head)
      pntr2=pointer+mhead-memloc(tmp_b)

* Determine the normalised coordinates of the point (0->NX-1, 0->NY-1)

      xj=(x-xmin)/dx
      yj=(y-ymin)/dy
      jx=int(xj)
      jy=int(yj)
      xj=xj-jx
      yj=yj-jy

* Determine coordinates of the lower left corner of the spline window

      jx=jx-nh
      jy=jy-nh

* Check if spline window is completely within grid

      if (jx.lt.-1 .or. jx+nwin.ge.nx .or. jy.lt.-1 .or. jy+nwin.ge.ny)
     |		then
         gridbspl=nan
	 return
      endif

* Load 6-by-6 subgrid

      if (ntype.eq.1) then
	 pntr2=pntr2+jy*nx+jx
         do ky=1,nwin
	    pntr2=pntr2+nx
            do kx=1,nwin
	       z(kx,ky)=tmp_b(pntr2+kx)
	    enddo
         enddo
      else if (ntype.eq.3) then
	 pntr2=pntr2/2+jy*nx+jx
         do ky=1,nwin
	    pntr2=pntr2+nx
            do kx=1,nwin
	       z(kx,ky)=tmp_s(pntr2+kx)
	    enddo
         enddo
      else if (ntype.eq.4) then
	 pntr2=pntr2/4+jy*nx+jx
         do ky=1,nwin
	    pntr2=pntr2+nx
            do kx=1,nwin
	       z(kx,ky)=tmp_i(pntr2+kx)
	    enddo
         enddo
      else if (ntype.eq.5) then
	 pntr2=pntr2/4+jy*nx+jx
         do ky=1,nwin
	    pntr2=pntr2+nx
            do kx=1,nwin
	       z(kx,ky)=tmp_f(pntr2+kx)
	    enddo
         enddo
      else
	 pntr2=pntr2/8+jy*nx+jx
         do ky=1,nwin
	    pntr2=pntr2+nx
            do kx=1,nwin
	       z(kx,ky)=tmp_d(pntr2+kx)
	    enddo
         enddo
      endif

* Conduct spline interpolation in horizontal direction

      do ky=1,nwin
         call gridbsplt(z(1,ky),nwin,w,u)
	 zy(ky)=gridbsplu(z(nh,ky),z(nh+1,ky),w(nh),w(nh+1),xj)
      enddo

* Conduct spline interpolation in vertical direction

      call gridbsplt(zy,nwin,w,u)
      gridbspl=gridbsplu(zy(nh),zy(nh+1),w(nh),w(nh+1),yj)*dz+z0
      end

      subroutine gridbsplt(y,n,w,u)
! Based on "spline" from Numerical Recipes
      integer*4 n,k
      real*8	y(n),w(n),u(n),p
      w(1)=0
      u(1)=0
      do k=2,n-1
         p=w(k-1)/2+2
	 w(k)=-0.5d0/p
	 u(k)=(3*(y(k+1)-2*y(k)+y(k-1))-u(k-1)/2)/p
      enddo
      w(n-1)=u(n-1)
      do k=n-2,n/2,-1
         w(k)=w(k)*w(k+1)+u(k)
      enddo
      end

      function gridbsplu(y0,y1,w0,w1,b)
! Based on "splint" from Numerical Recipes
      real*8	y0,y1,w0,w1,b,gridbsplu
      gridbsplu=y0+b*(y1-y0-w0/3-w1/6+b*(w0/2+b*(w1-w0)/6))
      end
