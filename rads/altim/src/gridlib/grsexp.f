      subroutine grsexp(a,b,nx,ny,nx2,ny2,ma,mb,error)
***********************************************************************
* Expand grid a into grid b. The grid is blown up from nx*ny to       *
* size nx2*ny2.                                                       *
*                                                                     *
* Variables:                                                          *
*            Input:                                                   *
*   a     (R4)  original array containing grid                        *
*   nx,ny (I4)  x- and y-dim of grid a                                *
* nx2,ny2 (I4)  x- and y-dim of grid b                                *
*   ma,mb (I4)  x-dimension of arrays a and b                         *
*                                                                     *
*            Output:                                                  *
*   b     (R4)  array containing object grid                          *
*   error (L4)  .true. on error                                       *
*                                                                     *
* Remko Scharroo, Delft, January 19, 1990.                            *
* With special thanks to Rene Zandbergen.                             *
***********************************************************************
      integer nx,ny,nx2,ny2,ma,mb,kx,ky
      real a(ma,*),b(mb,*)
      real xmul,ymul,fx,fy
      logical error
*
  550 format (a)
      error=.false.
*
      if (nx.lt.2 .or. nx.gt.ma .or. ny.lt.2 .or.
     . nx2.lt.2 .or. nx.gt.ma .or. ny2.lt.2) then
        write (0,550) 'grsexp: invalid array ranges received'
        goto 1399
      endif
      xmul=(nx2-1.)/(nx-1.)
      ymul=(ny2-1.)/(ny-1.)
      if (nx2.gt.mb) then
        write (0,550) 'grsexp: expanded dimension too large'
        goto 1399
      endif
      do ky=1,ny2
        fy=1.+(ky-1)/ymul
        do kx=1,nx2
          fx=1.+(kx-1)/xmul
          call grintp(a,fx,fy,ma,b(kx,ky),error)
          if (error) goto 1300
        enddo
      enddo
      return
*
 1300 write (0,550) 'grsexp: error passed'
 1399 error=.true.
      end
