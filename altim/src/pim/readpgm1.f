      subroutine readpgm1(line,i)
      integer*1 line(80)
      integer i,n,m
      do n=1,80
         call GRRFIL(i,1,line(n))
	 if (line(n).eq.10) goto 1
      enddo
1     do m=n,80
	 line(m)=32
      enddo
      end

      subroutine readpgm2(buf,pixmap,nx,ny,ci)
      integer*1 buf(nx,ny)
      integer pixmap(nx,ny),nx,ny,ci(0:255)
      integer ix,iy,i
      do iy=1,ny
	 do ix=1,nx
	    i=buf(ix,iy)
	    if (i.lt.0) i=i+256
	    pixmap(ix,ny-iy+1)=ci(i)
	 enddo
      enddo
      end
