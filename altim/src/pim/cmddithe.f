      subroutine cmddithe
      include "pim.inc"

* Set dithering for colours (dither1, def=.false.) and shades (dither2,
* def=.true.

      dither1=.false.
      dither2=.true.
      if (popcmd('DITHE',argum)) then
	 dither1=.true.
	 read (argum,*,iostat=ios) dither1,dither2
      endif

      end
