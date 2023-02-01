#ifdef DOUBLE
#define Real REAL*8
#else
#define Real REAL*4
#define GR3INT GR3INT4
#endif
**GR3INT -- Interpolate cubically in a two-dimensional Real array
*+
      FUNCTION  GR3INT (GRID, X, Y, NX, NY, MX)
      INTEGER*4 NX, NY, MX
      Real    GR3INT, GRID(MX,*), X, Y
*
* This subroutine interpolates cubically in a two-dimensional regular grid.
* The dimension of the grid is NX by NY. The real coordinates of the point to be
* interpolated are (X,Y) where 1 <= X <= NX and 1 <= Y <= NY.
*
* Arguments:
* GRID    (input): Two dimension array containing grid.
* X,Y     (input): Real grid coordinates of the point to be interpolated.
* NX,NY   (input): X- and Y- dimension of the grid.
* MX      (input): First dimension of the array as defined in the calling
*                  routine. MX >= NX.
* GR3INT (output): Cubically interpolated value at point (X,Y).
*-
* 29-Oct-1990 - Created: Pieter Visser and Remko Scharroo, Delft
* 17-May-1991 - New manual
*-----------------------------------------------------------------------
      Real const(16),trx,try,val
      integer*4 kx,ky,i,j,k,l
      integer*4 coef(16,16)
      data coef/
     . 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  1,
     . 0,  0,  0, -6,  0,  0,  3,  0,  0,  0,  0,  6,  3, -3,  0, -3,
     . 0,  6,  0,  3,  0,  0, -6,  0,  0, -6,  0, -3, -3,  6,  0,  3,
     . 0,  0,  0,  2,  0,  0,  3,  0,  0,  0,  0, -2,  1, -3,  0, -1,
     . 0,  0,  0, -6,  0,  0,  0,  3,  0,  0,  0,  3,  6,  0, -3, -3,
     . 0,  0,  0, 36,  0,  0,-18,-18,  9,  0,  0,-18,-18,  9,  9,  9,
     . 0,-36,  0,-18, 18,  0, 36,  9,-18, 18,  0,  9, 18,-18, -9, -9,
     . 0,  0,  0,-12,  0,  0,-18,  6,  9,  0,  0,  6, -6,  9,  3,  3,
     . 0,  0,  6,  3,  0,  0,  0, -6,  0,  0, -6, -3, -3,  0,  6,  3,
     . 0,  0,-36,-18,  0, 18,  9, 36,-18,  0, 18, 18,  9, -9,-18, -9,
     .36, 18, 18,  9,-36,-36,-18,-18, 36,-18,-18, -9, -9, 18, 18,  9,
     . 0,  0, 12,  6,  0, 18,  9,-12,-18,  0,  6, -6,  3, -9, -6, -3,
     . 0,  0,  0,  2,  0,  0,  0,  3,  0,  0,  0,  1, -2,  0, -3, -1,
     . 0,  0,  0,-12,  0,  0,  6,-18,  9,  0,  0, -6,  6,  3,  9,  3,
     . 0, 12,  0,  6, 18,  0,-12,  9,-18,  6,  0,  3, -6, -6, -9, -3,
     . 0,  0,  0,  4,  0,  0,  6,  6,  9,  0,  0,  2,  2,  3,  3,  1/
*
      kx=max(3,min(nx-1,int(x+1)))
      ky=max(3,min(ny-1,int(y+1)))
      trx=x-kx
      try=y-ky
*
      do i=1,16
         const(i)=0
      enddo
*
      l=0
      do i=-2,1
         do j=-2,1
            l=l+1
            val=grid(kx+i,ky+j)
            do k=1,16
               const(k)=const(k)+coef(k,l)*val
            enddo
         enddo
      enddo
*
      GR3INT=(const(1)+const(2)*trx+const(3)*try+const(4)*trx*try
     .       +const(5)*trx**2+const(6)*try**2
     .       +const(7)*trx*try**2+const(8)*try*trx**2
     .       +const(9)*trx**2*try**2
     .       +const(10)*trx**3+const(11)*try**3
     .       +const(12)*trx**3*try+const(13)*try*trx**3
     .       +const(14)*trx**3*try**2+const(15)*try**3*trx**2
     .       +const(16)*trx**3*try**3)/36
      end
