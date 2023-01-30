      subroutine cmdsymbo
      include "pim.inc"
      real s,dmin
      logical l,xy,pimopt,chropt,time

* Plot symbols at requested position

      call pgsave
360   if (pop1cmd('SYMBO',argum)) then
	 ci=c_fg
	 ls=1
	 lw=1
	 ch=def_ch
	 s=-1
	 dmin=0
	 xy=.true.
	 time=.false.
	 if (pimopt('-t',argum,dum,dum,dum,dum)) time=.true.
	 if (pimopt('-xy',argum,dum,dum,dum,dum)) xy=.true.
	 if (pimopt('-yx',argum,dum,dum,dum,dum)) xy=.false.
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('ch=',argum,ch,dum,dum,dum)
	 l=pimopt('symb=',argum,s,dum,dum,dum)
	 l=pimopt('mrk=',argum,s,dum,dum,dum)
	 l=pimopt('dmin=',argum,dmin,dum,dum,dum)
	 call pgsch(ch)
	 call pgsci(nint(ci))
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 write (0,'(a)') 'Plotting symbols ...'
	 if (chropt('file=',argum,name)) then
	    open (20,file=name,status='old')
	    call symbols(work1,maxgrd/2,nint(s),time,xy,dmin)
	    close (20)
	 else
	    read (argum,*,iostat=ios) x,y
	    call pmconv(1,x,y)
	    call pgpt(1,x,y,nint(s))
	 endif
	 goto 360
      endif
      call pgunsa
      end

      subroutine symbols(xy,n,sym,timecol,xyinput,dmin)
      integer n,i,sym
      logical xyinput,timecol
      real xy(n,2),xlast,ylast,dmin2,dmin,t
      character*80 rec
      include "pim.inc"
      i=0
      ios=0
      xlast=1e30
      ylast=1e30
      dmin2=dmin**2

1     continue
      read (20,'(a)',iostat=ios) rec
      if (ios.ne.0) then
	 call symbol1(i,xy(1,1),xy(1,2),sym)
	 return
      else if (rec(:1).eq.'#') then
         goto 1
      else if (xyinput) then
	 if (timecol) then
	    read (rec,*,iostat=ios) t,x,y
	 else
	    read (rec,*,iostat=ios) x,y
	 endif
      else
	 if (timecol) then
	    read (rec,*,iostat=ios) t,y,x
	 else
	    read (rec,*,iostat=ios) y,x
	 endif
      endif

      if ((x-xlast)**2+(y-ylast)**2.lt.dmin2) then
	 goto 1
      else if (x.lt.xw0 .or. x.gt.xw1 .or. y.lt.yw0 .or. y.gt.yw1) then
	 goto 1
      else if (i.eq.n) then
	 call symbol1(i,xy(1,1),xy(1,2),sym)
	 i=1
      else
	 i=i+1
      endif
      xy(i,1)=x
      xy(i,2)=y
      xlast=x
      ylast=y
      goto 1
      end


      subroutine symbol1(n,x,y,sym)
      integer n,sym
      real x(*),y(*)
      call pmconv(n,x,y)
      call pgpt(n,x,y,sym)
      end
