      subroutine cmdconto
      include "pim.inc"
      include "pimcont.inc"
      logical pixels,contours,labels,l,pimopt,cmdrange
      integer n

* CONTOUR

      call pgsave
      l=cmdrange()
200   if (pop1cmd('CONTO',cargum)) then
	 ci=c_cont
	 ls=1
	 lw=1
	 ch=def_ch
	 rintc=0
	 l=pimopt('ls=',cargum,ls,dum,dum,dum)
	 l=pimopt('lw=',cargum,lw,dum,dum,dum)
	 l=pimopt('ch=',cargum,ch,dum,dum,dum)
	 l=pimopt('ci=',cargum,ci,dum,dum,dum)
	 pixels=(inter.eq.2 .and. ls.eq.1 .and.
     &		lw.eq.1 .and. .not.dumped)
	 call strip1(cargum,cname)
 	 if (cname.eq.' '.or.cname.eq.'=') cname=sname
 	 if (cname.eq.sname .and. (lname.eq.sname .or. lname.eq.' '))
     &		then
 	    write (0,550) 'Copying contour grid from surface grid ...'
* This is a dummy action. WORK1 already contains surface grid and WORK2
* the interpolated version
 	 else
	    write (0,550) 'Reading contour grid ...'
	    call gridread(work1,work2,work4,ncx,ncy,cname,
     &		xc0,xc1,yc0,yc1)
	    if (pixels) then
	       write (0,550) 'Interpolating contour grid ...'
	       call pixinter(work1,ncx,ncy,xc0,xc1,yc0,yc1,work2,px,py)
	    endif
	 endif
	 fact=2.5
	 l=pimopt('sig=',cargum,fact,dum,dum,dum)
	 call grstat4(work1,ncx,ncy,ncx,yc0,yc1,fact,
     |		n,rminc,rmaxc,rmeanc,dum,rsigc)
	 call roundoff(rminc,rmaxc)
	 pgcint=max(ncx,ncy)/3
	 l=pimopt('rng=',cargum,rminc,rmaxc,rintc,dum)
	 l=pimopt('min=',cargum,rminc,dum,dum,dum)
	 l=pimopt('max=',cargum,rmaxc,dum,dum,dum)
	 l=pimopt('int=',cargum,rintc,dum,dum,dum)
	 l=pimopt('ci=',cargum,ci,dum,dum,dum)
	 labels=(pimopt('-l',cargum,dum,dum,dum,dum) .or.
     |		 pimopt('lab=',cargum,pgcint,dum,dum,dum))
	 contours=.not.pimopt('-n',cargum,dum,dum,dum,dum)
	 read (cargum,*,iostat=ios) rminc,rmaxc,rintc
	 if (rintc.eq.0) then
	    rintc=(rmaxs-rmins)/(nc1+1)
	    rminc=nint((rminc-rmins)/rintc-0.5)*rintc+rmins
	    rmaxc=nint((rmaxc-rmins)/rintc+0.5)*rintc+rmins
	 endif
	 write (0,'("... Selected range: ",2f12.6)') rminc,rmaxc

* If interpolation and a pixel map exists, and this is a pixel oriented device,
* do fast contouring

	 if (contours.and.pixels) then
	    write (0,550) 'Doing fast contouring ...'
	    call pgpxcont(px,py,px,work2,work3,work4,
     &		rminc,rmaxc,rintc,nint(ci))
	    contours=.false.
	 endif
	 if (.not.dumped) call dumppix
	 call pgsci(nint(ci))
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsch(ch)
	 call pimcont(work1,contours,labels)
	 goto 200
      endif
      call pgunsa
550   format (a)
      end
