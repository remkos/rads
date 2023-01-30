      subroutine boxes
      include 'xrgb.inc'
      integer i
      do i=0,maxind
	 call pgsci(i)
	 call pgsfs(1)
	 call pgrect(xl(i),xr(i),yb(i),yt(i))
	 call pgsci(1)
	 call pgsfs(2)
	 call pgrect(xl(i),xr(i),yb(i),yt(i))
      enddo
      end
