      subroutine grest1(a,b,fact,nx,ny,mx)
***********************************************************************
* Estimate factor to be applied to grid b in order to fit grid a.     *
*                                                                     *
* Variables:                                                          *
*            Input:                                                   *
*   a,b   (R4)  arrays with source grids                              *
*   nx,ny (I4)  x- and y-dim of grids                                 *
*   mx    (I4)  x-dimension of arrays a and b                         *
*                                                                     *
*            Output:                                                  *
*   fact  (R4)  multiplier to be applied to grid b                    *
*                                                                     *
* Remko Scharroo, Delft, January 24, 1990.                            *
* With special thanks to Rene Zandbergen.                             *
***********************************************************************
      integer kx,ky,nx,ny,mx
      real a(mx,*),b(mx,*)
      real sumab,sumbb,va,vb,fact
*
      sumab=0.
      sumbb=0.
      do ky=1,ny
         do kx=1,nx
            va=a(kx,ky)
            vb=b(kx,ky)
	    if (va+vb.lt.1e20) then
               sumab=sumab+va*vb
               sumbb=sumbb+vb*vb
	    endif
	 enddo
      enddo
      if (sumbb.eq.0.) then
         fact=0.
      else
         fact=sumab/sumbb
      endif
      end
