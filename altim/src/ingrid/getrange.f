      subroutine getrange(i0,i1,x0,x1,xmin,xmax,nx,quit,error)
      logical quit,error
      integer i0,i1,i,nx
      real xmin,xmax,x0,x1,fx/1.0/
      character*80 c

      error=.false.
      quit=.false.
550   format (a)

      read (5,550,end=190,err=20) c
      read (c,*,iostat=i) i0,i1
      if (i.ne.0) then
	 fx=(xmax-xmin)/(nx-1)
	 read (c,*,err=20) x0,x1
	 i0=nint(1+(x0-xmin)/fx)
	 i1=nint(1+(x1-xmin)/fx)
      endif
      if (i0.lt.1 .or. i0.gt.i1 .or. i1.gt.nx) then
	 write (0,550) 'ingrid: invalid range'
	 goto 20
      endif
      x0=xmin+fx*(i0-1)
      x1=xmin+fx*(i1-1)

      return

20    error=.true.
      return

190   quit=.true.
      return
      end


