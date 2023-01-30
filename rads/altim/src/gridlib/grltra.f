      subroutine grltra(a,bias,tilt,nx,ny,mx)
***********************************************************************
* Linear transformation of data in one grid (tilt'n'bias)             *
*                                                                     *
* Variables:                                                          *
*            Input/output:                                            *
*   a     (R4)  array with source/object grid                         *
*                                                                     *
*            Input:                                                   *
*   nx,ny (I4)  x- and y-dim of grid                                  *
*   tilt  (R4)  tilt to be multiplied with all points                 *
*   bias  (R4)  bias to be added to all points                        *
*   mx    (I4)  x-dimension of array a                                *
*                                                                     *
* Remko Scharroo, Delft, January 11, 1990.                            *
* With special thanks to Rene Zandbergen.                             *
***********************************************************************
      integer nx,ny,mx,kx,ky
      real a(mx,*)
      real tilt,bias
*
      do ky=1,ny
        do kx=1,nx
          a(kx,ky)=bias+tilt*a(kx,ky)
        enddo
      enddo
      end
