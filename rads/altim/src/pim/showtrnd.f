      subroutine showtrnd
      include "pim.inc"
      character*1 chr
      integer i,kx,ky
      real c(10)

* Use cursor to show trend in grids or XGF data

      if (.not.(gridtrnd.or.xgftrnd)) return

      call pgsave
      ch=def_ch
      ci=c_fg
      ls=1
      lw=1
      call pimopt('ch=',argum,ch,dum,dum,dum)
      call pimopt('ci=',argum,ci,dum,dum,dum)
      call pimopt('ls=',argum,ls,dum,dum,dum)
      call pimopt('lw=',argum,lw,dum,dum,dum)
      write (0,550) 'Click the left mouse button, or Q to quit ...'
30    call pgcurs(x,y,chr)
      if (chr.eq.'A') then
	 call pmcinv(1,x,y)
	 if (gridtrnd) then
	 call gridread(work1,work2,work4,nsx,nsy,sname,
     |          xs0,xs1,ys0,ys1)
	    kx=nint((x-xs0)/(xs1-xs0)*(nsx-1))+1
	    ky=nint((y-ys0)/(ys1-ys0)*(nsy-1))+1
	    do i=1,nstack
	       ty(i)=work2(kx+(ky-1)*nsx+(i-1)*nsx*nsy)
	    enddo
	    x=xs0+(kx-1)*(xs1-xs0)/(nsx-1)
	    y=ys0+(ky-1)*(ys1-ys0)/(nsy-1)
	    call gridstck(1,nstack,tx,ty,c)
	 else if (xgftrnd) then
	    call xgfstck(xname,x,y,nstack,tx,ty,c)
	 endif
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgtrend(x,y,c)
	 call pmconv(x,y)
	 goto 30
      endif
      call pgunsa
550   format (a)
      end
