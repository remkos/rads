**GRIDDY -- Derivative of grid in X derection
*+
      SUBROUTINE GRIDDY (H, G, NX, NY, MX)
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
      integer kx,ky,iy

* Combine derivative in y-direction

      do kx=1,nx
         g(kx,1)=(-3*h(kx,1)+4*h(kx,2)-h(kx,3))/2
         do ky=2,ny-1
            g(kx,ky)=(h(kx,ky+1)-h(kx,ky-1))/2
	 enddo
         g(kx,ny)=(3*h(kx,ny)-4*h(kx,ny-1)+h(kx,ny-2))/2
      enddo

* Mask out invalid fields

      do ky=1,ny
	 do kx=1,nx
	    if (abs(h(kx,ky)).ge.1e20) then
	       iy=max(2,min(ny-1,ky))
	       g(kx,iy-1)=1e30
	       g(kx,iy  )=1e30
	       g(kx,iy+1)=1e30
	    endif
	 enddo
      enddo
      end
