      subroutine cmdtitle
      include "pim.inc"

* Print title and adjust viewport

      call pgsave
      if (popcmd('TITLE',argum)) then
	 write (0,'(a)') 'Plotting title ...'
	 ci=c_fg
	 ls=1
	 lw=2
	 ch=1
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 call pimopt('ls=',argum,ls,dum,dum,dum)
	 call pimopt('lw=',argum,lw,dum,dum,dum)
	 call pimopt('ch=',argum,ch,dum,dum,dum)
	 call pgsci(nint(ci))
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsch(ch)
	 call pgmtext('B',1.5,0.5,0.5,argum)
	 yv0=yv0+0.05
	 call pgsvp(xv0,xv1,yv0,yv1)
      endif
      call pgunsa
      end
