*************************************************************
      subroutine rot90(ipos,a,b,mx,my)
*
* Rotate array a (mx,my) clockwise (ipos=1) or counter-clockwise (ipos=2)
* over 90 degrees. Object is array b (my,mx).
*
      integer mx,my,ipos,kx,ky
      real*4 a(mx,my),b(my,mx)
*
      if (ipos.eq.1) then
         do 100 ky=1,my
            do 100 kx=1,mx
  100          b(ky,mx-kx+1)=a(kx,ky)
      else
         do 110 ky=1,my
            do 110 kx=1,mx
  110          b(my-ky+1,kx)=a(kx,ky)
      endif
      end
