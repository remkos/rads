      subroutine initgrid
      implicit none
      include "grid.inc"
      integer grdidx

* Initialise xover rms grid statistics

      do grdidx=1,fnx*fny
         trms(grdidx)=0e0
         tnr(grdidx)=0
      enddo
      end

      subroutine wrgrid
      include "grid.inc"
      include "init.inc"
      integer ix,iy,det,ildet,undet
      integer px,bx,by,kx,ky,jx,jy,k,grdidx
      real*4  y,cx,w,value,weight
      real    frms(maxgrd)

* Start smoothing the grid

* px is used if the grid in wrapped over the +/-180 degree boundary,
* it is the number of grid cells covering 360 degrees.

      px=nint(360./fx)
      by=smcells

      grdidx=0
      det=0
      ildet=0
      undet=0
      do iy=1,fny
	 y=fy0+fy*(iy-1)
	 cx=cos(y*rad)
	 bx=nint(by/cx)
         do ix=1,fnx
	    grdidx=grdidx+1
	    value=0.
	    weight=0.

* We spread the information in latitude direction over 3 cells (-by...by),
* and in longitude direction over a strip (-bx...bx), where bx depends
* on the latitude (bx=2 for lat=48.2 to 66.4; bx=3 for lat=66.4 to 73.4; etc)
* The RMS contributions of each original cell are weighted by the distance
* to the cell that is evaluated. The weighting formula is:
* w = exp (-d**2) where d is the normalized distance with
* d**2 = (dx/fx/cos(lat))**2+(dy/fy)**2 = (kx*cx)**2 + ky**2

	    do ky=-by,by
	       jy=iy+ky
	       if (jy.ge.1 .and. jy.le.fny) then
	          do kx=-bx,bx
		     jx=mod(ix+kx-1+px,px)+1
		     if (jx.le.fnx) then
			k=(jy-1)*fnx+jx
		        w=exp(-(kx*cx)**2-ky**2)
		        value=value+w*trms(k)
		        weight=weight+w*tnr(k)
		     endif
		  enddo
	       endif
	    enddo

* If there are two few points in the neighbourhood, then weight is low,
* and an undetermined value is used. Otherwise, the sea surface variability
* is limited to the range 1-200 cm.

	    if (weight.ge.2.0) then
	       det=det+1
	       frms(grdidx)=sqrt(value/weight/2.)
	    else if (weight.ge.0.5) then
	       ildet=ildet+1
	       frms(grdidx)=sqrt(value/weight/2.)
	    else
	       undet=undet+1
	       frms(grdidx)=1e30
	    endif
	 enddo
      enddo

* Then write the grid

      open (10,file=res)
      close (10,status='delete')
      call gridwr4(res,fnx,fny,frms,fnx,fx0,fx1,fy0,fy1)

      end
