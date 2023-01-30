**ILLUMPIX -- Apply illumination grid to pixel map
*+
      SUBROUTINE ILLUMPIX (ILLUM, PIXMAP, N)
      INTEGER N
      REAL ILLUM(N)
      INTEGER PIXMAP(N)
*
* This routine recomputes the colour indices belonging to each pixel
* in the pixel map PIXMAP, depending on its current value (determined
* in SURF2PIX), and the illumination level in ILLUM. Also does the
* dithering in the shading direction.
*-
      real rest,vfact,v,cmax,c
      integer i,ic

      include "pim.inc"

      vfact=(nc2+1)/(rmaxl-rminl)
      rest=0.0
      cmax=nc2+0.999

      do i=1,n
	 ic=pixmap(i)
	 if (ic.eq.c_bad) then
	 else if (abs(illum(i)).gt.1e20) then
            pixmap(i)=c_bad
	 else
	    v=(illum(i)-rminl)*vfact
	    c=max(0.,min(cmax,v+rest))
	    ic=int(c)
	    pixmap(i)=sh_mat(pixmap(i),ic)
	    if (dither2) rest=c-ic+0.5
         endif
      enddo
      end
