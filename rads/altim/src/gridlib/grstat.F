#ifdef DOUBLE
#define Real REAL*8
#define GRSTAT GRSTAT8
#else
#define Real REAL*4
#define GRSTAT GRSTAT4
#endif
**GRSTAT -- Determine grid statistics
*+
      SUBROUTINE GRSTAT (GRID, NX, NY, MX, YMIN, YMAX, FACT,
     |           NVALID, ZMIN, ZMAX, ZMEAN, ZRMS, ZSIGMA)
      INTEGER*4 NX, NY, MX, NVALID
      Real GRID(MX,*), YMIN, YMAX, FACT,
     |		ZMIN, ZMAX, ZMEAN, ZRMS, ZSIGMA
*
* Subroutine to compute statistics of values contained on a grid.
* Weighting by latitude is invoked by specifying YMIN and YMAX as the
* minimum and maximum latitude of the grid. Weighting is disabled when
* both are set to zero.
*
* The routine also allow an iterative editing of outliers, which are off
* by more than FACT times the standard deviation from the mean. Use
* a non-positive value of FACT to disable editing.
*
* Upon return, the values ZMIN, ZMAX, ZMEAN, ZSIGMA are set to the minimum,
* maximum, mean, and standard deviation of the NVALID values in the grid
* (after editing).
*
* Input arguments:
*  GRID      : Grid of dimension NX*NY. The first dimension of the array
*              as defined in the calling (sub)program is MX.
*  NX, NY    : Number of grid values in X and Y direction.
*  MX        : Dimension of array GRID as defined in the calling (sub)program.
*  YMIN,YMAX : Latitude boudaries of the grid in degrees. To disable
*              weighting by the cosine of latitude, specify 0 for both values.
*  FACT      : Factor determines the editing range:
*              [ ZMEAN - FACT*ZSIGMA ; ZMEAN + FACT*ZSIGMA ]
*
* Output arguments:
*  NVALID    : Number of non-edited grid values.
*  ZMIN,ZMAX : Minimum and maximum value in the grid.
*  ZMEAN     : (Weighted) mean of the grid values.
*  ZRMS      : (Weighted) RMS of the grid values.
*  ZSIGMA    : (Weighted) standard deviation of the grid values.
*-
* 07-Nov-1996 - New version
*-----------------------------------------------------------------------
      Real    z,lo,hi,w,dy,rad
      real*8  wtot0,wtot,sum1,sum2
      integer kx,ky,iter

      wtot=0
      iter=0
      lo=-1e20
      hi= 1e20
      rad=atan(1e0)/45
      dy=(ymax-ymin)/max(1,ny-1)

   10 wtot0=wtot
      wtot=0
      nvalid=0
      zmin=+1e30
      zmax=-1e30
      sum1=0
      sum2=0
      iter=iter+1
      do ky=1,ny
         w=cos((ymin+(ky-1)*dy)*rad)
         do kx=1,nx
            z=grid(kx,ky)
            if (z.ge.lo .and. z.le.hi) then
               nvalid=nvalid+1
               wtot=wtot+w
               sum1=sum1+z*w
               sum2=sum2+z*z*w
               zmin=min(zmin,z)
               zmax=max(zmax,z)
            endif
         enddo
      enddo
      zmean=sum1/wtot
      zrms=sqrt(sum2/wtot)
      zsigma=sqrt(sum2/wtot-(sum1/wtot)**2)

      lo=zmean-fact*zsigma
      hi=zmean+fact*zsigma

      if (fact.gt.0 .and. wtot.ne.wtot0 .and. iter.lt.5) goto 10

      end
