      subroutine addrgb(idx,i,j)
      integer idx,i,j
      include 'xrgb.inc'

      changed=.true.
      cmap(i,idx)=cmap(i,idx)+j
      if (cmap(i,idx).gt.255) cmap(i,idx)=255
      if (cmap(i,idx).lt.0) cmap(i,idx)=0
      end
