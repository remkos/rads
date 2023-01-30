      subroutine copy(n,a,b)
      integer n,i
      real a(n),b(n)

      do i=1,n
	 b(i)=a(i)
      enddo
      end
