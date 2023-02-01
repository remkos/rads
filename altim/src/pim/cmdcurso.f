      subroutine cmdcurso(rint)
      include "pim.inc"
      real z,dx,dy,gr1int4,rint
      integer unit
      character*1 chr

* Use cursor

      call pgsave
      if (popcmd('CURSO',argum)) then
         ch=def_ch
         ci=c_cont
         ls=1
         lw=1
         call pimopt('ch=',argum,ch,dum,dum,dum)
         call pimopt('ci=',argum,ci,dum,dum,dum)
         call pimopt('ls=',argum,ls,dum,dum,dum)
         call pimopt('lw=',argum,lw,dum,dum,dum)
         call pimopt('int=',argum,rint,dum,dum,dum)
	 write (0,550) 'Re-loading contour grid ...'
	 call gridread(work2,work3,work4,ncx,ncy,cname,
     &		xc0,xc1,yc0,yc1)
	 dx=(xc1-xc0)/(ncx-1)
	 dy=(yc1-yc0)/(ncy-1)
	 if (argum.eq.' ') then
	    unit=6
	 else
	    unit=20
	    open (unit,file=argum,status='unknown')
10	    read (unit,550,end=20) name
	    goto 10
20	    continue
	 endif
	 write (0,550) 'Click the left mouse button, or Q to quit ...'
30	 call pgcurs(x,y,chr)
40	 format (2f8.3,1x,f8.3)
	 if (chr.eq.'A') then
	    call pmcinv(1,x,y)
	    z=gr1int4(work2,(x-xc0)/dx+1,(y-yc0)/dy+1,ncx,ncy,ncx)
	    write (unit,40) x,y,z
	    call pgsls(nint(ls))
	    call pgslw(nint(lw))
	    call annotate(x,y,work2,ncx,ncy,rint,nint(ci))
	    goto 30
	 endif
	 if (unit.ne.6) close (unit)
      endif
      call pgunsa
550   format (a)
      end
