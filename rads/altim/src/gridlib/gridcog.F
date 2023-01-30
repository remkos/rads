#ifdef DOUBLE
#define Real REAL*8
#define GRIDCOGn GRIDCOG8
#else
#define Real REAL*4
#define GRIDCOGn GRIDCOG4
#endif
**GRIDCOGn -- Find the Centre-Of-Gravity in a Real grid
*+
      SUBROUTINE GRIDCOGn (Z, NX, NY, MX, ZMIN, ZMAX, XCOG, YCOG)
      INTEGER*4 NX, NY, MX
      Real    Z(MX,NY), ZMIN, ZMAX, XCOG(NY+1), YCOG(NX+1)

* This routine termines the Centre-Of-Gravity (COG) of a two dimensional
* matrix (grid) of Real values. Along with it, it determined the COG along
* each gridline in both directions.
*
* The COG along a vertical gridline I is defined as SUM Z(I,J)*J / SUM Z(I,J),
* where J runs from 1 to NY. A similar definition holds for the COG along a
* horizontal gridline. Values of Z beyond the limits ZMIN and ZMAX are not
* taken into account.
*
* The position of the COG in vertical direction is stored in YCOG,
* horizontal in XCOG.
*
* The total COG of the matrix is found by taking into account all
* gridlines. The X and Y coordinates are stored as XCOG(NY+1) and
* YCOG(NX+1).
*
* Coordinates of the COG are expressed in gridpoints, counting from 1 to NX
* or 1 to NY for horizontal and vertical coordinates respectively. A value
* of 0 is stored when the COG could not be computed (division by zero).
*
* Arguments:
*  Z     (input): Input matrix (grid)
*  NX    (input): Number of gridpoints in X (horizontal) direction
*  NY    (input): Number of gridpoints in Y (vertical) direction
*  MX    (input): First dimension of Z as determined in the calling
*                 program or routine
*  ZMIN  (input): Lower limit of Z for editing
*  ZMAX  (input): Upper limit of Z for editing
*  XCOG (output): Horizontal positions of the COG. XCOG(J) contains the
*                 value for vertical gridline J, XCOG(NY+1) for the entire
*                 matrix
*  YCOG (output): Vertical positions of the COG. YCOG(I) contains the
*                 value for horizontal gridline I, YCOG(NX+1) for the entire
*                 matrix
*-
* 17-Apr-1999 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      real*8 prod,sum,t_prod,t_sum
      integer*4 kx,ky

* Compute horizontal COG positions

      t_prod=0d0
      t_sum =0d0

      do ky=1,ny
         prod=0d0
         sum=0d0
         do kx=1,nx
            if (z(kx,ky).ge.zmin .and. z(kx,ky).le.zmax) then
               prod=prod+z(kx,ky)*kx
               sum =sum +z(kx,ky)
            endif
         enddo
         if (sum.eq.0d0) then
            xcog(ky)=0d0
         else
            xcog(ky)=prod/sum
         endif
         t_prod=t_prod+prod
         t_sum =t_sum +sum
      enddo
      if (t_sum.eq.0d0) then
         xcog(ny+1)=0d0
      else
         xcog(ny+1)=t_prod/t_sum
      endif

* Compute vertical COG positions

      t_prod=0d0
      t_sum =0d0

      do kx=1,nx
         prod=0d0
         sum=0d0
         do ky=1,ny
            if (z(kx,ky).ge.zmin .and. z(kx,ky).le.zmax) then
               prod=prod+z(kx,ky)*ky
               sum =sum +z(kx,ky)
            endif
         enddo
         if (sum.eq.0d0) then
            ycog(kx)=0d0
         else
            ycog(kx)=prod/sum
         endif
         t_prod=t_prod+prod
         t_sum =t_sum +sum
      enddo
      if (t_sum.eq.0d0) then
         ycog(nx+1)=0d0
      else
         ycog(nx+1)=t_prod/t_sum
      endif

      return
      end
