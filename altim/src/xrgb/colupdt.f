      subroutine colupdt(i)
      integer i
      include 'xrgb.inc'

      if (cmap(1,i).ge.0) then
         call pgscr(i,cmap(1,i)/255.,cmap(2,i)/255.,cmap(3,i)/255.)
	 if (.not.used(i)) then
	    call pgsch(0.75)
	    call pgsci(1)
	    call pgptext
     .		(xl(i)+.5*xblksize,yb(i)+.1*yblksize,0.0,0.5,ctxt(i))
	    call pgsch(1.00)
	 endif
	 used(i)=.true.
      else
         call pgscr(i,cmap(1,0)/255.,cmap(2,0)/255.,cmap(3,0)/255.)
	 if (used(i)) then
	    call pgsch(0.75)
	    call pgsci(0)
	    call pgptext
     .		(xl(i)+.5*xblksize,yb(i)+.1*yblksize,0.0,0.5,ctxt(i))
	    call pgsch(1.00)
	 endif
	 used(i)=.false.
      endif

      end
