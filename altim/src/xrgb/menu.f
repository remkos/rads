      subroutine menu(i,text1)
      integer i
      character*(*) text1
      include 'xrgb.inc'
      if (i.eq.-1) then
         call pgsci(0)
         call pgsfs(1)
         call pgrect(0.,1.,0.,1.)
         call pgsci(1)
         call pgsfs(2)
         call pgsch(1.0)
      else if (i.eq.0) then
         call pgsci(0)
         call pgsfs(1)
         call pgrect(0.,1.,0.,.245)
         call pgsci(1)
         call pgsfs(2)
         call pgsch(1.0)
      else if (i.ge.13) then
	 call pgsci(0)
	 call pgsfs(1)
         call pgrect(0.,1.,0.,0.07)
	 call pgsci(1)
	 call pgsfs(2)
	 call pgptext(0.5,0.03,0.0,0.5,text1)
      else
	 call pgrect(bx0(i),bx1(i),by0(i),by1(i))
	 call pgptext(bx0(i)+dbx/2,by0(i)+dby/3,0.0,0.5,text1)
      endif
      end
