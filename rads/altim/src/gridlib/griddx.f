**GRIDDX -- Derivative of grid in X derection
*+
      SUBROUTINE GRIDDX (H, G, NX, NY, MX)
      INTEGER NX, NY, MX
      REAL H(MX,*), G(MX,*)
*
* Compute derivative of geometry grid H along the X axis
*
* Arguments:
*  H     (input): Geometry grid
*  G    (output): Derivative grid (DH/DX)
*  NX,NY (input): Dimension of the grid
*  MX    (input): First dimension of arrays H and G
*-
*  3-Jan-1995 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer kx,ky,ix

* Compute derivative in x-direction

      do ky=1,ny
         g(1,ky)=(-3*h(1,ky)+4*h(2,ky)-h(3,ky))/2
         do kx=2,nx-1
            g(kx,ky)=(h(kx+1,ky)-h(kx-1,ky))/2
	 enddo
         g(nx,ky)=(3*h(nx,ky)-4*h(nx-1,ky)+h(nx-2,ky))/2
      enddo

* Mask out invalid fields

      do ky=1,ny
	 do kx=1,nx
	    if (abs(h(kx,ky)).ge.1e20) then
	       ix=max(2,min(nx-1,kx))
	       g(ix-1,ky)=1e30
	       g(ix  ,ky)=1e30
	       g(ix+1,ky)=1e30
	    endif
	 enddo
      enddo
      end
