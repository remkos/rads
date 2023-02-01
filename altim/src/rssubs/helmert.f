      real*8 xyzfrom(3,100),xyzto(3,100),sigma(100),coeff(7),
     |	xyzto2(3,100),rms(2)
      integer i,j,n
      i=1
      read (*,*) coeff
10    read (*,*,end=20) (xyzfrom(j,i),j=1,3)
      read (*,*) sigma(i)
      i=i+1
      goto 10
20    n=i-1
      write (*,*) coeff
      call helmert2(coeff,n,xyzfrom,xyzto)
      call helmert1(n,xyzfrom,xyzto,sigma,coeff,rms)
      call helmert2(coeff,n,xyzfrom,xyzto2)
      write (*,*) coeff,rms
      do i=1,n
         write (*,*) (xyzto2(j,i)-xyzto(j,i),j=1,3)
      enddo
      end
