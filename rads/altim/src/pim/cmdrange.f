      function cmdrange()
      logical cmdrange
      include "pim.inc"

* Overrule or set (in the absence of grids) the height range

      if (popcmd('RANGE',argum)) then
	 read (argum,*,iostat=ios) rmins,rmaxs
	 cmdrange=.true.
      else
	 cmdrange=.false.
      endif

      end
