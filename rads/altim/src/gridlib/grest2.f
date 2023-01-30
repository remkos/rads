      subroutine grest2(a,b,fb,c,fc,nx,ny,mx,sing)
***********************************************************************
* Estimate linear combination of grids b and c, in order to fit grid a*
*                                                                     *
* Variables:                                                          *
*            Input:                                                   *
*   a,b,c (R4)  arrays with source grids                              *
*   nx,ny (I4)  x- and y-dim of grids                                 *
*   mx    (I4)  x-dimension of arrays a, b and c                      *
*                                                                     *
*            Output:                                                  *
*   sing  (L4)  .true. if b and c are linear combinations             *
*   fa,fb (R4)  multipliers to be applied to grid a and b             *
*                                                                     *
* Remko Scharroo, Delft, January 24, 1990.                            *
* With special thanks to Rene Zandbergen.                             *
***********************************************************************
      integer nx,ny,mx,kx,ky
      real a(mx,*),b(mx,*),c(mx,*)
      real sumab,sumac,sumbb,sumbc,sumcc,va,vb,fb,fc,vc,codet
      logical sing
*
      sumbb=0.
      sumbc=0.
      sumcc=0.
      sumab=0.
      sumac=0.
*
      do ky=1,ny
         do kx=1,nx
            va=a(kx,ky)
            vb=b(kx,ky)
            vc=c(kx,ky)
	    if (va+vb+vc.lt.1e20) then
               sumbb=sumbb+vb*vb
               sumbc=sumbc+vb*vc
               sumcc=sumcc+vc*vc
               sumab=sumab+va*vb
               sumac=sumac+va*vc
	    endif
         enddo
      enddo
*
      codet=sumbb*sumcc-sumbc*sumbc
      sing=(codet.eq.0.0)
      if (sing) then
         fb=0.0
         fc=0.0
      else
         fb=(sumab*sumcc-sumac*sumbc)/codet
         fc=(sumbb*sumac-sumbc*sumab)/codet
      endif
      end
