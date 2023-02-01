      subroutine untag(i)
* Untag the selected colour
      integer i
      include 'xrgb.inc'
      call pgslw(3)
      call pgsci(i)
      call pgsfs(2)
      call pgrect(xl(i)+.2*xblksize,xr(i)-.2*xblksize,
     .            yb(i)+.5*yblksize,yt(i)-.2*yblksize)
      call pgslw(1)
      call pgsci(1)
      call pgsfs(1)
      end
