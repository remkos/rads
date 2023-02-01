**SURF2PIX -- Convert surface grid to pixel map
*+
      SUBROUTINE SURF2PIX (SURF, PIXMAP, N)
      INTEGER N
      REAL SURF(N)
      INTEGER PIXMAP(N)
*
* This routine assignes colour indices to the surface levels in the
* array SURF. The output is a pixel map PIXMAP. Also dithering is
* performed when requested.
*-
      real rest,vfact,v,cmax,c
      integer i,ic

      include "pim.inc"

      vfact=(nc1+1)/(rmaxs-rmins)
      rest=0.0
      cmax=nc1+0.999

      do i=1,n
	 c=surf(i)
	 if (.not.slog) then
	 else if (c.le.0 .or. c.gt.1e20) then
	    c=1e35
	 else
	    c=log10(c)
	 endif
         if (abs(c).gt.1e20) then
            pixmap(i)=c_bad
	 else
            v=(c-rmins)*vfact
            c=c_0+max(0.,min(cmax,v+rest))
            ic=int(c)
            pixmap(i)=ic
            if (dither1) rest=c-ic-0.5
         endif
      enddo
      end
