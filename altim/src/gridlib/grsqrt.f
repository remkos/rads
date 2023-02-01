      subroutine grsqrt(a,nx,ny,mx)
***********************************************************************
* Take sqaure root of data in one grid                                *
*                                                                     *
* Variables:                                                          *
*            Input/output:                                            *
*   a     (R4)  array with source/object grid                         *
*                                                                     *
*            Input:                                                   *
*   nx,ny (I4)  x- and y-dim of grid                                  *
*   mx    (I4)  x-dimension of array a                                *
*                                                                     *
* Remko Scharroo, Delft, January 11, 1990.                            *
* With special thanks to Rene Zandbergen.                             *
***********************************************************************
      integer kx,ky,nx,ny,mx
      real a(mx,*)
*
      do ky=1,ny
        do kx=1,nx
          a(kx,ky)=sqrt(abs(a(kx,ky)))
        enddo
      enddo
      end
