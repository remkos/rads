*************************************************************
      subroutine subarr(b,nxb,nyb,c,nxc,mx,my)
*
* Substitute array b (nxb,nyb) in array c (nxc,*). The lower left
* corner of b (1,1) will be placed at position (mx,my) in array c.
*
      integer nxb,nyb,nxc,mx,my,kx,ky
      real*4 b(nxb,*),c(nxc,*)

      do 100 kx=1,nxb
         do 100 ky=1,nyb
  100       c(kx+mx-1,ky+my-1)=b(kx,ky)
      end
