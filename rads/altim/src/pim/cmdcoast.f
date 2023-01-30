      subroutine cmdcoast(layer)
      integer layer
      include "pim.inc"   

* Plot coast lines

      call pgsave
  320 if (pop1cmd('COAST',argum)) then
	 ci=c_coast
	 ls=1
	 lw=1
	 x=1
         name='1.cil'
	 call pimopt('layer=',argum,x,dum,dum,dum)
	 if (layer.ne.nint(x)) goto 320
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 call pimopt('ls=',argum,ls,dum,dum,dum)
	 call pimopt('lw=',argum,lw,dum,dum,dum)
	 call strip(argum,name)
	 read (argum,*,iostat=ios) lw
	 call pgsci(nint(ci))
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 write (0,'(a)') 'Drawing coastlines ...'
	 call pmwdb(name,1,0)
	 goto 320
      endif
      call pgunsa
      end
