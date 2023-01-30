      integer buflen
      parameter (buflen=200000)
      character*18 xgfbuf(buflen)
      integer rec0,rec1,fd
      common /xgfcom/ xgfbuf,fd,rec0,rec1

      integer maxgrd,nx,ny
      parameter (maxgrd=2000000)
      real*8 x0,x1,dx,y0,y1,dy,h(maxgrd)
      common /cdiswgrid/ x0,x1,dx,y0,y1,dy,h,nx,ny
