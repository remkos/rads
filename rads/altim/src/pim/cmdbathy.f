      subroutine cmdbathy
      include "pim.inc"
      real bath(4)
      integer i

* Plot bathymetry

      call pgsave
  300 if (pop1cmd('BATHY',argum)) then
	 write (0,'(a)') 'Drawing bathymetry ...'
	 ls=1
	 lw=1
	 ci=c_cont
         bath(1)=-2000
         bath(2)=-3500
         bath(3)=-1
         bath(4)=-1
         name='2.nat.bth'
	 call pimopt('ls=',argum,ls,dum,dum,dum)
	 call pimopt('lw=',argum,lw,dum,dum,dum)
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 call pimopt('lev=',argum,bath(1),bath(2),bath(3),bath(4))
	 call strip(argum,name)
	 read (argum,*,iostat=ios) lw,bath
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsci(nint(ci))
	 do i=1,4
	    if (nint(bath(i)).ne.-1) then
	       call pmwdb(name,nint(bath(i)),nint(bath(i)))
	    endif
	 enddo
	 goto 300
      endif
      call pgunsa
      end
