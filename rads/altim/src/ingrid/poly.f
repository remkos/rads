*************************************************************
      subroutine poly(x,z,n,pol)
*
* Compute polynomial of x into z.
*
      include "ingrid.inc"

      real pol(MPOL),x,z
      integer n,k,inv

      z=0
      if (n.lt.0 .or. n.ge.MPOL) return
      do 100 k=0,n
         inv=n+1-k
  100    z=x*z+pol(inv)
      end
