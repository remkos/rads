      subroutine cmdrect
      include "pim.inc"   
      real xx(2),yy(2)
      integer i,n

* Plot a pie chunk

      call pgsave
  320 if (pop1cmd('RECT',argum)) then
	 ci=c_fg
	 fs=1
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 call pimopt('fs=',argum,fs,dum,dum,dum)
	 read (argum,*,iostat=ios) xx,yy
	 call pmconv(2,xx,yy)
	 call pgsci(nint(ci))
	 call pgsfs(nint(fs))
	 write (0,'(a)') 'Drawing rectangle ...'
	 call pgrect(xx(1),xx(2),yy(1),yy(2))
	 goto 320
      endif
      call pgunsa
      end
