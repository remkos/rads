      integer n,i
      real*8 x, y(10), g(10), yp

      read(*,*) n,x
      read(*,*) (y(i),i=1,2*n)
      call everett(n,x,y,yp,g)
      write (*,*) (g(i),i=1,2*n)
      write (*,*) yp
      end

      SUBROUTINE EVERETT(N, X, Y, YP, G)
      INTEGER*4 N, M
      REAL*8	X, Y(1-N:N), YP, G(2,0:N-1)

* Arguments:
*  N   (input): Order of the interpolater = half the number of data points
*  X   (input): Point at which the interpolated value is required
*  Y   (input): Table of input values Y(1-N), Y(2-N), ... Y(N)
*  YP (output): Interpolated value
*  G  (output): Central differences

      integer*4 i,j,k,r,r2,bin
      real*8	s

      g(1,0)=y(0)
      g(2,0)=y(1)
      do r=1,n-1
	 k=0
	 bin=1
	 g(1,r)=y(-r)  +y(r)
	 g(2,r)=y(-r+1)+y(r+1)
	 r2=2*r+1
         do j=-r+1,r-1
	    k=k+1
	    bin=-bin*(r2-k)/k
	    write(*,*) r,j,bin
	    g(1,r)=g(1,r)+bin*y(j)
	    g(2,r)=g(2,r)+bin*y(j+1)
	 enddo
      enddo

      yp=0
      s=1
      k=0
      do r=0,n-1
         yp=yp+s*g(1,r)
	 k=k+1
	 s=s*x/k
	 yp=yp+s*g(2,r)
	 k=k+1
	 s=s*(1-x)/k
         write (*,*) yp
      enddo

      end
