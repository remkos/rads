#ifdef DOUBLE
#define Real REAL*8
#else
#define Real REAL*4
#define GR2INT GR2INT4
#endif
**GR2INT -- Interpolate quadratically in a two-dimensional array
*+
      FUNCTION  GR2INT (GRID, X, Y, NX, NY, MX)
      INTEGER*4 NX, NY, MX
      Real    GR2INT, GRID(MX,*), X, Y
*
* This subroutine interpolates quadratically in a two-dimensional regular grid.
* The dimension of the grid is NX by NY. The real coordinates of the point to be
* interpolated are (X,Y) where 1 <= X <= NX and 1 <= Y <= NY.
*
* Arguments:
* GRID    (input): Two dimension array containing grid.
* X,Y     (input): Real grid coordinates of the point to be interpolated.
* NX,NY   (input): X- and Y- dimension of the grid.
* MX      (input): First dimension of the array as defined in the calling
*                  routine. MX >= NX.
* GR2INT (output): Quadratically interpolated value at point (X,Y).
*-
* 29-Oct-1990 - Created: Pieter Visser and Remko Scharroo, Delft
* 17-May-1991 - New manual
*  9-Dec-1994 - Common creation of GR2INT and GR2INT4
*-----------------------------------------------------------------------
      Real const(9),trx,try,val
      integer*4 kx,ky,i,j,k,l
      integer*4 coef(9,9)
      data coef/
     .0, 0, 0, 1, 0, 0,-1,-1, 1,
     .0,-2, 0, 0, 2, 0, 2, 0,-2,
     .0, 0, 0,-1, 0, 0,-1, 1, 1,
     .0, 0,-2, 0, 0, 2, 0, 2,-2,
     .4, 0, 0, 0,-4,-4, 0, 0, 4,
     .0, 0, 2, 0, 0, 2, 0,-2,-2,
     .0, 0, 0,-1, 0, 0, 1,-1, 1,
     .0, 2, 0, 0, 2, 0,-2, 0,-2,
     .0, 0, 0, 1, 0, 0, 1, 1, 1/
*
      kx=max(2,min(nx-1,nint(x)))
      ky=max(2,min(ny-1,nint(y)))
      trx=x-kx
      try=y-ky
*
      do i=1,9
         const(i)=0
      enddo
*
      l=0
      do i=-1,1
         do j=-1,1
            l=l+1
            val=grid(kx+i,ky+j)
            do k=1,9
               const(k)=const(k)+coef(k,l)*val
            enddo
         enddo
      enddo
*
      GR2INT=(const(1)+const(2)*trx+const(3)*try+const(4)*trx*try
     .       +const(5)*trx**2+const(6)*try**2
     .       +const(7)*trx*try**2+const(8)*try*trx**2
     .       +const(9)*trx**2*try**2)/4
      end
