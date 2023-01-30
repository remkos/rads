      subroutine cmdland(layer)
      integer layer
      include "pim.inc"

* Plot land mask

      call pgsave
  310 if (pop1cmd('LAND',argum)) then
	 ci=c_land
	 x=1
         name='1.lnd'
	 call pimopt('layer=',argum,x,dum,dum,dum)
	 if (layer.ne.nint(x)) goto 310
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 if (argum.ne.' ') name=argum
	 write (0,'(a)') 'Drawing land mask ...'
	 call pmwdb(name,nint(ci),c_bad)
	 goto 310
      endif
      call pgunsa
      end
