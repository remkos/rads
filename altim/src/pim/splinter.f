**SPLINTER -- Interpolate grid to higher resolution
*
      subroutine splinter(a,nx,ny,b,px,py,w,nw)

      integer nx,ny,px,py,nw
      real*4 a(nx,ny),b(px,py),w(nw)
      integer naux,i,k
      real f

      k=0
      do i=1,nx
         k=k+1
         w(k)=i-1
      enddo
      do i=1,ny
         k=k+1
         w(k)=i-1
      enddo
      f=real(nx-1)/real(px-1)
      do i=1,px
         k=k+1
         w(k)=(i-1)*f
      enddo
      f=real(ny-1)/real(py-1)
      do i=1,py
         k=k+1
         w(k)=(i-1)*f
      enddo
      naux=nw-nx-ny-px-py

      call scsin2(w(1),w(nx+1),a,nx,ny,nx,w(nx+ny+1),w(nx+ny+px+1),
     &            px,py,b,px,w(nx+ny+px+py+1),naux)
      return
      end
