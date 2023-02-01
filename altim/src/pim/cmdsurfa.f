      subroutine cmdsurfa
      include "pim.inc"
      logical l,pimopt,cmdrange
      integer n

* SURFACE

      if (sname.ne.' ') then
	 write (0,550) 'Reading surface grid ...'
	 call gridread(work1,work2,work4,nsx,nsy,sname,
     |		xs0,xs1,ys0,ys1)
	 fact=2.5
         l=pimopt('sig=',sargum,fact,dum,dum,dum)
	 call grstat4(work1,nsx,nsy,nsx,ys0,ys1,fact,
     |		n,rmins,rmaxs,rmeans,dum,rsigs)
	 write (0,'("... Grid statistics:",5f12.6)')
     |		rmins,rmaxs,rmeans,sqrt(rmeans**2+rsigs**2),rsigs
	 call gridcog4(work1,nsx,nsy,nsx,-1e20,1e20,xscog,yscog)
	 write (0,'("... Grid centre-of-gravity:",2f12.6)')
     |		xs0+(xscog(nsy+1)-1)*(xs1-xs0)/(nsx-1),
     |		ys0+(yscog(nsx+1)-1)*(ys1-ys0)/(nsy-1)
	 call roundoff(rmins,rmaxs)
	 l=cmdrange()
	 slog=pimopt('-log',sargum,dum,dum,dum,dum)
	 l=pimopt('rng=',sargum,rmins,rmaxs,dum,dum)
	 l=pimopt('min=',sargum,rmins,dum,dum,dum)
	 l=pimopt('max=',sargum,rmaxs,dum,dum,dum)
	 gridtrnd=pimopt('-show',sargum,dum,dum,dum,dum)
	 read (sargum,*,iostat=ios) rmins,rmaxs
         if (l) then
	    rminl=rmeans-fact*rsigs
	    rmaxl=rmeans+fact*rsigs
	    call roundoff(rminl,rmaxl)
	    rmins=max(rminl,rmins)
	    rmaxs=min(rmaxl,rmaxs)
         endif
	 write (0,'("... Selected range: ",2f12.6)') rmins,rmaxs
	 if (inter.ge.1) then
            write (0,550) 'Interpolating surface grid ...'
	    call pixinter(work1,nsx,nsy,xs0,xs1,ys0,ys1,work2,px,py)
	    write (0,550) 'Computing colour indices ...'
	    call surf2pix(work2,work3,px*py)
	 else
	    write (0,550) 'Computing colour indices ...'
	    call surf2pix(work1,work3,px*py)
	 endif
	 dumped=.false.
      endif
550   format (a)
      end
