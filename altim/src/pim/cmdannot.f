      subroutine cmdannot(rint)
      include "pim.inc"
      real rint

* Annotate contours

      call pgsave
      if (popcmd('ANNOT',argum)) then
	 ch=def_ch
         lw=1
	 ls=1
         ci=c_cont
	 call pimopt('ls=',argum,ls,dum,dum,dum)
	 call pimopt('ch=',argum,ch,dum,dum,dum)
	 call pimopt('ci=',argum,ci,dum,dum,dum)
	 call pimopt('lw=',argum,lw,dum,dum,dum)
	 call pimopt('int=',argum,rint,dum,dum,dum)
	 write (0,550) 'Annotating contours ...'
	 open (20,file=argum,form='formatted',status='old',err=999)
	 call pgsch(ch)
         call pgslw(nint(lw))
         call pgsls(nint(ls))
	 rewind (20)
  270       read (20,550,end=275) argum
            if (argum(1:1).eq.'#') goto 270
	    read (argum,*) x,y
	    call annotate(x,y,work2,px,py,rint,nint(ci))
	    goto 270
  275    continue
	 close (20)
      endif
  550 format (a)
  999 call pgunsa
      end
