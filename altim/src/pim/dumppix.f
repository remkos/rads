      subroutine dumppix
      include "pim.inc"
      real xblc,xtrc,yblc,ytrc

* Dump the pixel map.
* If INTER is on, the pixels run throughout the entire PLOTAREA.
* If INTER is off, the pixels overlay the GRIDAREA boundaries by half.
* If INTER is prevented, do a boundary coordinate transformation first.

      write (0,"('Plotting pixel map ...')")
      if (inter.ge.1) then
         call pgqwin(xblc,xtrc,yblc,ytrc)
      else
	 call pmconv(1,xs0,ys0)
	 call pmconv(1,xs1,ys1)
	 xblc=xs0-(xs1-xs0)/(px-1)/2
	 xtrc=xs1+(xs1-xs0)/(px-1)/2
	 yblc=ys0-(ys1-ys0)/(py-1)/2
	 ytrc=ys1+(ys1-ys0)/(py-1)/2
      endif
      call pgpixl(work3,px,py,1,px,1,py,xblc,xtrc,yblc,ytrc)
      dumped=.true.
      end
