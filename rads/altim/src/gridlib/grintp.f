      subroutine grintp(z,fx,fy,mx,val,error)
***********************************************************************
* Interpolate value in point (fx,fy) somewhere in 2-dimensional grid  *
*                                                                     *
* Variables:                                                          *
*            Input:                                                   *
*   z     (R4)  array containing grid                                 *
*   fx,fy (R4)  real coordinates of point within grid                 *
*               1-eps <= fx <= mx+eps   and   1-eps <= fy             *
*   mx    (I4)  x-dimension of array z                                *
*                                                                     *
*            Output:                                                  *
*   val   (R4)  linearly interpolated value in point (fx,fy)          *
*   error (L4)  .true. if nx or ny out of bounds                      *
*                                                                     *
* Remko Scharroo, Delft, January 19, 1990.                            *
* With special thanks to Rene Zandbergen.                             *
***********************************************************************
      integer kx,ky,mx
      real*4 z(mx,*),fx,fy,frx1,frx2,fry1,fry2,eps,val
      logical error
*
  550 format (a)
      error=.false.
      eps=1e-5
*
      if (fx.lt.1-eps .or. fx.gt.mx+eps) then
        write (0,550) 'grintp: illegal x request'
        goto 1301
      endif
      if (fy.lt.1-eps) then
        write (0,550) 'grintp: illegal y request'
        goto 1301
      endif
      kx=max(1,int(fx-eps))
      ky=max(1,int(fy-eps))
      frx2=fx-real(kx)
      fry2=fy-real(ky)
      frx1=1.-frx2
      fry1=1.-fry2
      val=fry1*(frx1*z(kx,ky)+frx2*z(kx+1,ky))
      val=val+fry2*(frx1*z(kx,ky+1)+frx2*z(kx+1,ky+1))
      return
 1301 error=.true.
      end
