**INTAB2 -- Interpolate a table using 2nd order polynomial (parabola)
*+
      SUBROUTINE INTAB2 (NDIM, NVEC, TABLE, T1, TN, T, VEC)
      INTEGER*4 NDIM, NVEC
      REAL*8  TABLE(NDIM,NVEC), T1, TN, T, VEC(NDIM)
*
* This subroutine interpolates a NDIM dimensional vector VEC at a given
* time T from a table of NVEC (>=3) vectors, equally spaced in time, starting
* at T1 and ending at TN. 'Time' can be seen in this respect as any
* independent variable of which the vector is a function. A 2nd order
* polynomial interpolation is used to interpolate the vector.
*
* Arguments:
*  NDIM  (input) : Dimension (= number of elements) of the vector to be
*                  interpolated. Thus NDIM=1 denotes the interpolation of a
*                  scalar.
*  NVEC  (input) : Number of table entries. The first table entry is for
*                  time T1, the last for time TN. NVEC must be greater or
*                  equal to 3.
*  TABLE (input) : Array containing NVEC vectors of dimension NDIM.
*  T1    (input) : Time corresponding to the first vector stored in TABLE
*                  ( TABLE(1,1)...TABLE(NDIM,1) )
*  TN    (input) : Time corresponding to the last vector stored in TABLE
*                  ( TABLE(1,NVEC)...TABLE(NDIM,NVEC) )
*  T     (input) : Time at which a vector must be interpolated.
*  VEC  (output) : Interpolated NDIM-dimensional vector at time T. The dimension
*                  defined in the calling (sub)program must be at least NDIM.
*-
* 10-Mar-1992 - Created [Remko Scharroo]
*-----------------------------------------------------------------------
      REAL*8  TREL
      INTEGER*4 ITREL
      TREL=(T-T1)/(TN-T1)*(NVEC-1)+1
      ITREL=MAX(1,MIN(NINT(TREL)-1,NVEC-2))
      CALL INTER2(NDIM,TREL-ITREL,TABLE(1,ITREL),VEC)
      END

      SUBROUTINE INTER2(NDIM,T,X,XT)
      INTEGER*4 NDIM, I
      REAL*8 T,X(NDIM,3),XT(NDIM),DT,DT2,A,B,C
      DT=T-1
      DT2=DT*DT
      DO 10 I=1,NDIM
	 A=X(I,2)
	 B=(X(I,3)-X(I,1))/2
	 C=X(I,3)-A-B
	 XT(I)=A+B*DT+C*DT2
   10 CONTINUE
      END
