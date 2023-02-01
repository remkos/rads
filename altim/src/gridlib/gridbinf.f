**GRIDBINF -- Return information about buffered grid
*+
      SUBROUTINE GRIDBINF (POINTER, NX, NY, XMIN, XMAX, YMIN, YMAX,
     |                     ZMIN, ZMAX)
      INTEGER*4 POINTER, NX, NY
      REAL*8    XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX

* This routine returns information about a grid previously loaded with
* the subroutine GRIDBUFF.
*
* Input argument:
*  POINTER    : Pointer to grid structure as returned by GRIDBUFF
*
* Output arguments:
*  NX, NY     : Number of gridpoints in X- and Y-direction
*  XMIN, XMAX : X-coordinate of the left- and rightmost gridpoint
*  YMIN, YMAX : Y-coordinate of the lower- and uppermost gridpoint
*  ZMIN, ZMAX : Minimum and maximum value in the grid
*-
* 26-Jul-2006 - Avoiding use of %val
* 22-Jul-2005 - Use allocated memory
* 21-Nov-2000 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4 mhead
      parameter (mhead=104)
      integer*4	ihdr(mhead/4)
      real*8	dhdr(mhead/8)
      equivalence (ihdr,dhdr)

      call memget(pointer,mhead,ihdr)
      nx  =ihdr(3)
      ny  =ihdr(4)
      xmin=dhdr(5)
      xmax=dhdr(6)
      ymin=dhdr(8)
      ymax=dhdr(9)
      zmin=dhdr(11)
      zmax=dhdr(12)
      end
