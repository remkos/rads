      subroutine cmdillum
      include "pim.inc"
      real scale,angle
      logical l,l2,pimopt
      integer n

* ILLUMINATION

      if (lname.ne.' ') then
 	 if (lname.eq.sname) then
 	    write (0,550)
     .		'Copying illumination grid from surface grid ...'
* This is a dummy action. WORK1 already contains
* surface grid
 	 else
	    write (0,550) 'Reading illumination grid ...'
	    call gridread(work1,work2,work4,nlx,nly,lname,
     &		xl0,xl1,yl0,yl1)
	 endif

	 angle=135
	 shade0=0.20
	 shade1=0.95

	 l=pimopt('ang=',largum,angle,dum,dum,dum)
	 write (0,550) 'Computing slopes ...'
*	 call grillu(work1,work4,angle,nlx,nly,nlx)
*	 scale=(xl1-xl0)/(nlx-1)
*	 do ix=1,nlx*nly
*	    work4(ix)=work4(ix)/scale
*	 enddo
	 call gridillu(work1,work5,work4,work4,angle)
	 fact=2.5
         l2=pimopt('sig=',largum,fact,dum,dum,dum)
	 call grstat4(work4,nlx,nly,nlx,yl0,yl1,fact,
     |		n,rminl,rmaxl,rmeanl,dum,rsigl)
	 write (0,'("... Grid statistics:",5f12.6)')
     |          rminl,rmaxl,rmeanl,sqrt(rmeanl**2+rsigl**2),rsigl
	 l=pimopt('rng=',largum,rminl,rmaxl,dum,dum)
	 l=pimopt('min=',largum,rminl,dum,dum,dum)
	 l=pimopt('max=',largum,rmaxl,dum,dum,dum)
	 l=pimopt('shade=',largum,shade0,shade1,dum,dum)
	 read (largum,*,iostat=ios) rminl,rmaxl
	 if (l2) then
            rminl=max(rmeanl-fact*rsigl,rminl)
            rmaxl=min(rmeanl+fact*rsigl,rmaxl)
         endif
	 write (0,'("... Selected range: ",2f12.6)') rminl,rmaxl
	 if (inter.ge.1) then
            write (0,550)'Interpolating slopes ...'
	    call pixinter(work4,nlx,nly,xl0,xl1,yl0,yl1,work5,px,py)
	    write (0,550) 'Recomputing colour indices ...'
	    call illumpix(work5,work3,px*py)
	 else
	    write (0,550) 'Recomputing colour indices ...'
	    call illumpix(work4,work3,px*py)
	 endif
	 dumped=.false.
      endif
550   format (a)
      end
