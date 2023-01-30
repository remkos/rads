**GRILLU -- 'Illuminate' geometry grid
*+
      SUBROUTINE GRILLU (H, G, ANGLE, NX, NY, MX)
      INTEGER NX, NY, MX
      REAL H(MX,*), G(MX,*), ANGLE
*
* Illuminate a geometry grid H in direction ANGLE and obtain grid G
*
* Arguments:
*  H     (input): Geometry grid
*  G    (output): Illumination grid = derivative of H along a line
*                 which makes an angle ANGLE with the 'east-west' axis
*  ANGLE (input): Angle of illumination (in degrees), 0 is in 'east-west'
*                 direction, 90 is in a 'north-south' direction
*  NX,NY (input): Dimension of the grid
*  MX    (input): First dimension of arrays H and G
*-
* 19-Jan-1990 - Created by Remko Scharroo from version by Rene Zandbergen
*  3-Jan-1995 - New manual
*-----------------------------------------------------------------------
      real s,c,z
      integer kx,ky,ix,iy

      z=angle*atan(1.)/45
      c=cos(z)/2
      s=sin(z)/2

* Compute derivative in x-direction (multiplied by factor 2)

      do ky=1,ny
         g(1,ky)=3*h(1,ky)-4*h(2,ky)+h(3,ky)
         do kx=2,nx-1
            g(kx,ky)=h(kx-1,ky)-h(kx+1,ky)
	 enddo
         g(nx,ky)=-3*h(nx,ky)+4*h(nx-1,ky)-h(nx-2,ky)
      enddo

* Combine derivative in y-direction with above computed derivative in x-dir.

      do kx=1,nx
         g(kx,1)=c*g(kx,1)+s*(3*h(kx,1)-4*h(kx,2)+h(kx,3))
         do ky=2,ny-1
            g(kx,ky)=c*g(kx,ky)+s*(h(kx,ky-1)-h(kx,ky+1))
	 enddo
         g(kx,ny)=c*g(kx,ny)+
     .             s*(-3*h(kx,ny)+4*h(kx,ny-1)-h(kx,ny-2))
      enddo

* Mask out invalid fields

      do ky=1,ny
	 do kx=1,nx
	    if (abs(h(kx,ky)).ge.1e20) then
	       ix=max(2,min(nx-1,kx))
	       iy=max(2,min(ny-1,ky))
	       g(ix  ,iy  )=1e30
	       g(ix  ,iy-1)=1e30
	       g(ix+1,iy  )=1e30
	       g(ix  ,iy+1)=1e30
	       g(ix-1,iy  )=1e30
	    endif
	 enddo
      enddo
      end
