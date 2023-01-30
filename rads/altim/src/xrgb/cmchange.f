      subroutine cmchange(idx)
      include 'xrgb.inc'
      character*80 text
      logical log(3),hls/.false./
      integer idx,i,j,chls(3)
      real h,l,s

  552 format(a1,1x,i3.3)
      if (idx.lt.0) then
      call menu(0,' ')
      call menu(12,'Done')
      call menu(13,'Click colour to change')
1     call pgcurse(x,y,ch)
      call whatindx(idx)
      ich=ichar(ch)
      if (ich.ne.65) then
	 goto 1
      else if (idx.lt.0) then
	 return
      endif
      endif
c
c Make RGB menu
c
      call menu(0,' ')
      call menu(1,'+')
      call menu(5,'+')
      call menu(9,'+')
      call menu(3,'-')
      call menu(7,'-')
      call menu(11,'-')
      call menu(4,'HLS')
      call menu(8,'RGB')
      call menu(12,'Done')
      call menu(13,
     .'Click + or - with left bottom (slow) or right (fast)')

    2 log(1)=.true.
      log(2)=.true.
      log(3)=.true.
      call tag(idx)

    3 call pgscf(1)
      call pgsch(1.2)
      call grxhls(cmap(1,idx)/255.,cmap(2,idx)/255.,cmap(3,idx)/255.,
     .h,l,s)
      chls(1)=nint(h)
      chls(2)=nint(l*100)
      chls(3)=nint(s*100)
      do i=1,3
	    j=i*4-2
	    call pgsci(0)
	    call pgsfs(1)
	    call pgrect(bx0(j),bx1(j),by0(j),by1(j))
	    call pgsfs(2)
	    call pgsci(1)
	    if (hls) then
               write (text,552) texthls(i),chls(i)
	    else
               write (text,552) textrgb(i),cmap(i,idx)
	    endif
	    call menu(j,text)
            log(i)=.false.
      enddo
      call pgsch(1.)

    4 call pgcurse(x,y,ch)
      ich=ichar(ch)
      if (ich.ne.65 .and. ich.ne.88) goto 4
      call whatindx(j)
      j=-j
      i=j/4+1
      if (j.eq.4) then
	 hls=.true.
	 goto 3
      else if (j.eq.8) then
	 hls=.false.
	 goto 3
      else if (mod(j,4).eq.1) then
c Add 1 or 10
         if (hls) then
	    if (ich.eq.65) call addhls(idx,i,1)
	    if (ich.eq.88) call addhls(idx,i,10)
	 else
	    if (ich.eq.65) call addrgb(idx,i,1)
	    if (ich.eq.88) call addrgb(idx,i,10)
	 endif
         log(i)=.true.
	 call colupdt(idx)
	 goto 3
      else if (mod(j,4).eq.3) then
c Subtract 1 or 10
         if (hls) then
	    if (ich.eq.65) call addhls(idx,i,-1)
	    if (ich.eq.88) call addhls(idx,i,-10)
	 else
	    if (ich.eq.65) call addrgb(idx,i,-1)
	    if (ich.eq.88) call addrgb(idx,i,-10)
	 endif
         log(i)=.true.
	 call colupdt(idx)
	 goto 3
      endif
      call untag(idx)
      call whatindx(idx)
      if (idx.ge.0) goto 2
      end
