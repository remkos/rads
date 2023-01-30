      subroutine cmdline
      include "pim.inc"
      real xa,xb,ya,yb,dmin,offx,offy
      logical l,xy,pimopt,chropt

* Plot symbols at requested position

      call pgsave
360   if (pop1cmd('LINE',argum)) then
	 ci=c_fg
	 ls=1
	 lw=1
	 dmin=1e15
	 xy=.true.
	 offx=0
	 offy=0
	 if (pimopt('-xy',argum,dum,dum,dum,dum)) then
	    xy=.true.
	 endif
	 if (pimopt('-yx',argum,dum,dum,dum,dum)) then
	    xy=.false.
	 endif
	 l=pimopt('cut=',argum,dmin,dum,dum,dum)
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('offxy=',argum,offx,offy,dum,dum)
	 l=pimopt('offyx=',argum,offy,offx,dum,dum)
	 call pgsci(nint(ci))
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 write (0,'(a)') 'Plotting line ...'
	 if (chropt('file=',argum,name)) then
	    open (20,file=name,status='old')
	    call line(work1,maxgrd/2,xy,offx,offy,dmin)
	    close (20)
	 else
	    read (argum,*,iostat=ios) xa,ya,xb,yb
	    call pmconv(1,xa,ya)
	    call pmconv(1,xb,yb)
	    call pgmove(xa,ya)
	    call pgdraw(xb,yb)
	 endif
	 goto 360
      endif
      call pgunsa
      end

      subroutine line(xy,n,xyinput,offx,offy,dmin)
      integer n,i
      logical xyinput
      real xy(n,2),dmin2,dmin,offx,offy
      character*80 rec

      include "pim.inc"
      i=0
      ios=0
      dmin2=dmin**2

1     continue
      read (20,'(a)',iostat=ios) rec
      if (ios.ne.0) then
	 call line1(i,xy(1,1),xy(1,2))
	 return
      else if (rec.eq.' ') then
	 call line1(i,xy(1,1),xy(1,2))
	 goto 1
      else if (rec(:1).eq.'#') then
         goto 1
      else if (xyinput) then
	 read (rec,*,iostat=ios) x,y
      else
	 read (rec,*,iostat=ios) y,x
      endif
      x=x+offx
      y=y+offy

      if (i.eq.0) then
      else if (i.eq.n) then
	 call line1(i,xy(1,1),xy(1,2))
      else if ((x-xy(i,1))**2+(y-xy(i,2))**2.gt.dmin2) then
	 call line1(i,xy(1,1),xy(1,2))
      endif

      i=i+1
      xy(i,1)=x
      xy(i,2)=y
      goto 1
      end

      subroutine line1(n,x,y)
      integer n
      real x(*),y(*)
      call pmconv(n,x,y)
      call pgline(n,x,y)
      n=0
      end
