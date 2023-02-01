      function cmdbad()
      logical cmdbad
      include "pim.inc"

* Set lower and upper bound of values that are not set to invalid.
* Anything below the lower bound and above the upper bound it set to invalid.

      rminb=-1e20
      rmaxb=1e20
      if (popcmd('BAD',argum)) then
	 read (argum,*,iostat=ios) rminb,rmaxb
	 cmdbad=.true.
      endif

      end
