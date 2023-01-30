      subroutine cmdbox
      include "pim.inc"
      character*20 xopt,yopt
      real xint,yint
      integer ixsub,iysub

* Draw box

      call pgsave
      ls=1
      lw=1
      ci=c_fg
      ch=def_ch
      xopt='bcnsti'
      yopt='bcnstiv'
      xint=0
      yint=0
      ixsub=0
      iysub=0
      if (popcmd('BOX',argum)) then
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 call pimopt('ls=',argum,ls,dum,dum,dum)
	 call pimopt('lw=',argum,lw,dum,dum,dum)
	 call pimopt('ch=',argum,ch,dum,dum,dum)
	 read (argum,*,iostat=ios)
     .		lw,xopt,xint,ixsub,yopt,yint,iysub
      endif
      if (lw.eq.0) return
      write (0,'(a)') 'Plotting box ...'
      call pgsci(nint(ci))
      call pgsch(ch)
      call pgsls(nint(ls))
      call pgslw(nint(lw))
      call pmbox(xopt,xint,ixsub,yopt,yint,iysub)
      call pgunsa
      end
