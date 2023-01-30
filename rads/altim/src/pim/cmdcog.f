      subroutine cmdcog
      include "pim.inc"
      real s,dsx,dsy
      logical l,pimopt,xmark,ymark,omark
      integer kx,ky,sxmark,symark,somark,step

* Plot symbols for centre-of-gravity
* Options:
* -x: only x-positions (default symbol 2)
* -y: only y-positions (default symbol 7)
* -o: only total COG (default symbol 891)
* Further: lw= ch= ls= ci= symb= step=

      call pgsave
360   if (pop1cmd('COG',argum)) then
	 ci=c_fg
	 ls=1
	 lw=1
	 ch=def_ch
	 xmark=.true.
	 ymark=.true.
	 omark=.true.
	 sxmark=2
	 symark=7
	 somark=8
	 step=1
	 if (pimopt('-x',argum,dum,dum,dum,dum)) then
	    ymark=.false.
	    omark=.false.
	 endif
	 if (pimopt('-y',argum,dum,dum,dum,dum)) then
	    xmark=.false.
	    omark=.false.
	 endif
	 if (pimopt('-o',argum,dum,dum,dum,dum)) then
	    xmark=.false.
	    ymark=.false.
	 endif
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('ch=',argum,ch,dum,dum,dum)
	 if (pimopt('step=',argum,s,dum,dum,dum)) step=s
	 if (pimopt('symb=',argum,s,dum,dum,dum).or.
     |	     pimopt('mrk=',argum,s,dum,dum,dum)) then
	    sxmark=nint(s)
	    symark=nint(s)
	    somark=nint(s)
	 endif
	 call pgsch(ch)
	 call pgsci(nint(ci))
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 dsx=(xs1-xs0)/(nsx-1)
	 dsy=(ys1-ys0)/(nsy-1)
	 if (xmark) then
	    write (0,'(a)') 'Plotting COG symbols (x) ...'
	    do ky=1,nsy,step
	       if (xscog(ky).ge.1) then
	          x=xs0+dsx*(xscog(ky)-1)
	          y=ys0+dsy*(ky-1)
	          call pmconv(1,x,y)
	          call pgpt(1,x,y,sxmark)
	       endif
	    enddo
	 endif
	 if (ymark) then
	    write (0,'(a)') 'Plotting COG symbols (y) ...'
	    do kx=1,nsx,step
	       if (yscog(kx).ge.1) then
	          x=xs0+dsx*(kx-1)
	          y=ys0+dsy*(yscog(kx)-1)
	          call pmconv(1,x,y)
	          call pgpt(1,x,y,symark)
	       endif
	    enddo
	 endif
	 if (omark) then
	    write (0,'(a)') 'Plotting COG symbol  (o) ...'
	    if (xscog(nsy+1).ge.1 .and. yscog(nsx+1).ge.1) then
	       x=xs0+dsx*(xscog(nsy+1)-1)
	       y=ys0+dsy*(yscog(nsx+1)-1)
	       call pmconv(1,x,y)
	       call pgpt(1,x,y,somark)
	    endif
	 endif
	 goto 360
      endif
      call pgunsa
      end
