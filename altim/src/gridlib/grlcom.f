      subroutine grlcom(a,fa,b,fb,c,nx,ny,mx)
***********************************************************************
* Make a grid in c the linear combination of two grids in a and b.    *
*                                                                     *
* Variables:                                                          *
*            Input:                                                   *
*   a,b   (R4)  arrays with source grids                              *
*   fa,fb (R4)  multipliers to be applied to grid a and b             *
*   nx,ny (I4)  x- and y-dim of grids                                 *
*   mx    (I4)  x-dimension of arrays a and b                         *
*                                                                     *
*            Output:                                                  *
*   c     (R4)  array containing object grid                          *
*                                                                     *
* Remko Scharroo, Delft, January 19, 1990.                            *
* With special thanks to Rene Zandbergen.                             *
***********************************************************************
      integer kx,ky,nx,ny,mx
      real a(mx,*),b(mx,*),c(mx,*)
      real fa,fb
*
      do ky=1,ny
        do kx=1,nx
          c(kx,ky)=fa*a(kx,ky)+fb*b(kx,ky)
        enddo
      enddo
      end
