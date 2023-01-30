**INTAB8 -- Interpolate a table using 8th order Legendre interpolation
*+
      SUBROUTINE INTAB8 (NDIM, NVEC, TABLE, T1, TN, T, VEC)
      INTEGER*4 NDIM, NVEC
      REAL*8  TABLE(NDIM,NVEC), T1, TN, T, VEC(NDIM)
*
* This subroutine interpolates a NDIM dimensional vector VEC at a given
* time T from a table of NVEC (>=9) vectors, equally spaced in time, starting
* at T1 and ending at TN. 'Time' can be seen in this respect as any
* independent variable of which the vector is a function. A 8th order
* Legendre interpolation is used to interpolate the vector.
*
* Arguments:
*  NDIM  (input) : Dimension (= number of elements) of the vector to be
*                  interpolated. Thus NDIM=1 denotes the interpolation of a
*                  scalar.
*  NVEC  (input) : Number of table entries. The first table entry is for
*                  time T1, the last for time TN. NVEC must be greater or
*                  equal to 9.
*  TABLE (input) : Array containing NVEC vectors of dimension NDIM.
*  T1    (input) : Time corresponding to the first vector stored in TABLE
*                  ( TABLE(1,1)...TABLE(NDIM,1) )
*  TN    (input) : Time corresponding to the last vector stored in TABLE
*                  ( TABLE(1,NVEC)...TABLE(NDIM,NVEC) )
*  T     (input) : Time at which a vector must be interpolated.
*  VEC  (output) : Interpolated NDIM-dimensional vector at time T. The dimension
*                  defined in the calling (sub)program must be at least NDIM.
*-
* 20-Jan-1992 - Created: Remko Scharroo, DUT/SOM (c)
*-----------------------------------------------------------------------
      REAL*8  TREL
      INTEGER*4 ITREL
      TREL=(T-T1)/(TN-T1)*(NVEC-1)+1
      ITREL=MAX(1,MIN(NINT(TREL)-4,NVEC-8))
      CALL INTER8(NDIM,TREL-ITREL,TABLE(1,ITREL),VEC)
      END
