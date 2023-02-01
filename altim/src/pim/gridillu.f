      subroutine gridillu(z,zx,zy,zi,angle)
      include "pim.inc"

      real z(nlx,nly),zx(nlx,nly),zy(nlx,nly),zi(nlx,nly)
      real a(3),b(3),dx,dy,cosa,sina,rad,angle,lx,ly
      integer kx,ky

      call griddx(z,zx,nlx,nly,nlx)
      call griddy(z,zy,nlx,nly,nlx)

      dx=(xl1-xl0)/(nlx-1)
      dy=(yl1-yl0)/(nly-1)
      rad=atan(1e0)/45
      cosa=-cos(angle*rad)/dx
      sina=-sin(angle*rad)/dy
      if (project.eq.1) then
      do ky=1,nly
         do kx=1,nlx
	    zi(kx,ky)=zx(kx,ky)*cosa+zy(kx,ky)*sina
	 enddo
      enddo
      else
      do ky=1,nly
         do kx=1,nlx
	    a(1)=xl0+(kx-1)*dx
	    b(1)=yl0+(ky-1)*dy
	    a(2)=a(1)+dx
	    b(2)=b(1)
	    a(3)=a(1)
	    b(3)=b(1)+dy
	    call pmconv(3,a,b)
	    a(2)=a(2)-a(1)
	    b(2)=b(2)-b(1)
	    a(3)=a(3)-a(1)
	    b(3)=b(3)-b(1)
	    lx=sqrt(a(2)**2+b(2)**2)
	    ly=sqrt(a(3)**2+b(3)**2)
	    if (lx.gt.1e-6 .and. ly.gt.1e-6) then
	       zi(kx,ky)=
     .		zx(kx,ky)*(a(2)*cosa+b(2)*sina)/lx +
     .		zy(kx,ky)*(a(3)*cosa+b(3)*sina)/ly
	    else
	       zi(kx,ky)=1e30
	    endif
	 enddo
      enddo
      endif
      end
