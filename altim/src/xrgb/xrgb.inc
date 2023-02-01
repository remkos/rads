* xrgb.inc
*
      integer maxind
      parameter (maxind=144)
      character*1 ch
      integer ich
      real x,y
      common /pointer/ x,y,ich,ch
      real bx0(12),bx1(12),by0(12),by1(12),dbx,dby
      common /menuboxes/ bx0,bx1,by0,by1,dbx,dby
      real xl(0:maxind),xr(0:maxind),yb(0:maxind),yt(0:maxind)
      real xblksize,yblksize
      character*1 textrgb(3),texthls(3)
      common /colorboxes/ xl,xr,yb,yt,xblksize,yblksize,textrgb,
     .texthls
      integer*2 cmap(3,0:255)
      character*80 cmfile
      character*3 ctxt(0:255)
      logical used(0:maxind)
      common /color/ cmap,used,cmfile,ctxt
      logical changed
      common /parameters/ changed
