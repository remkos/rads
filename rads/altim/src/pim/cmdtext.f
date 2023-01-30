      subroutine cmdtext
      include "pim.inc"
      real angle,just

* Plot text at requested position

      call pgsave
350   if (pop1cmd('TEXT',argum)) then
	 ci=c_fg
	 ls=1
	 lw=1
	 ch=def_ch
	 angle=0
	 just=0.5
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 call pimopt('ls=',argum,ls,dum,dum,dum)
	 call pimopt('lw=',argum,lw,dum,dum,dum)
	 call pimopt('ch=',argum,ch,dum,dum,dum)
	 call pimopt('just=',argum,just,dum,dum,dum)
	 call pimopt('angle=',argum,angle,dum,dum,dum)
	 call pgsch(ch)
	 call pgsci(nint(ci))
	 call pgslw(nint(ls))
	 call pgslw(nint(lw))
	 write (0,'(a)') 'Plotting text ...'
	 if (argum(1:5).eq.'file=') then
	    open (20,file=argum(6:),status='old')
351	    read (20,'(a)',end=355) argum
	    read (argum,*,iostat=ios) x,y,name,angle,just
	    call pmconv(1,x,y)
	    call pgptxt(x,y,angle,just,name)
	    goto 351
355	    close (20)
	 else
	    read (argum,*,iostat=ios) x,y,name,angle,just
	    call pmconv(1,x,y)
	    call pgptxt(x,y,angle,just,name)
	 endif
	 goto 350
      endif
      call pgunsa
      end
