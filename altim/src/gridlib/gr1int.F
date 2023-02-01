#ifdef DOUBLE
#define Real REAL*8
#define BIG 1d30
#else
#define Real REAL*4
#define BIG 1e30
#define GR1INT GR1INT4
#define GR1IN1 GR1IN14
#endif
**GR1INT -- Interpolate linearly in a two-dimensional array
*+
      FUNCTION  GR1INT (Z, X, Y, NX, NY, MX)
      INTEGER*4 NX, NY, MX
      Real    GR1INT, Z(MX,*), X, Y
*
* This subroutine interpolates linearly in a two-dimensional regular grid.
* The dimension of the grid is NX by NY. The real coordinates of the point to be
* interpolated are (X,Y) where 1 <= X <= NX and 1 <= Y <= NY.
*
* Arguments:
* Z       (input): Two dimension array containing grid.
* X,Y     (input): Real grid coordinates of the point to be interpolated.
* NX,NY   (input): X- and Y- dimension of the grid.
* MX      (input): First dimension of the array as defined in the calling
*                  routine. MX >= NX.
* GR1INT (output): Linearly interpolated value at point (X,Y).
*-
* 22-May-1991 - Created - Remko Scharroo.
* 16-Nov-1994 - Triangular interpolation around invalid points
*  9-Dec-1994 - Common creation of GR1INT and GR1INT4
* 22-Nov-2000 - New interpolation scheme
*-----------------------------------------------------------------------
      Real GR1IN1
      integer*4 kx,ky,jx,jy
*
* Return 'invalid' value if requested location lies outside the grid
*
      GR1INT=BIG
      if (x.lt.1 .or. x.gt.nx .or.
     &    y.lt.1 .or. y.gt.ny) return

      kx=min(nx-1,int(x))
      ky=min(ny-1,int(y))
      jx=kx+1
      jy=ky+1

      GR1INT=GR1IN1(z(kx,ky),z(jx,ky),z(kx,jy),z(jx,jy),
     |	x-kx,y-ky)
      end

      function GR1IN1(zll,zlr,zul,zur,a,b)
      Real GR1IN1,a,b,wll,wlr,wul,wur,wtot,zll,zlr,zul,zur

* Determine the weights for all corner points. If undetermined, a
* corner point does not contribute to the interpolation.
* (a,b) is the relative position of (x,y) wrt (kx,ky)

      wll=0; wlr=0; wul=0; wur=0
      if (zll.lt.1d20) wll=(1-a)*(1-b)
      if (zlr.lt.1d20) wlr=   a *(1-b)
      if (zul.lt.1d20) wul=(1-a)* b
      if (zur.lt.1d20) wur=   a * b
      wtot=wll+wlr+wul+wur

* Total weight should at least be 0.5 (nomally 1) for interpolation

      if (wtot.gt.0.5d0) then
         GR1IN1=(zll*wll+zlr*wlr+zul*wul+zur*wur)/wtot
      else
         GR1IN1=BIG
      endif
      end
