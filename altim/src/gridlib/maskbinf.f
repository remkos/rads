**MASKBINF -- Return information about buffered mask
*+
      SUBROUTINE MASKBINF (POINTER, NX, NY, XMIN, XMAX, YMIN, YMAX)
      INTEGER*4 POINTER, NX, NY
      REAL*8    XMIN, XMAX, YMIN, YMAX

* This routine returns information about a mask previously loaded with
* the subroutine MASKBUFF.
*
* Input argument:
*  POINTER    : Pointer to grid structure as returned by GRIDBUFF
*
* Output arguments:
*  NX, NY     : Number of maskcells in X- and Y-direction
*  XMIN, XMAX : X-coordinate of the outer limits of
*               left- and rightmost maskcell
*  YMIN, YMAX : Y-coordinate of the outer limits of
*               lower- and uppermost maskcell
*-----------------------------------------------------------------------
      real*8 zmin,zmax
      call gridbinf (pointer,nx,ny,xmin,xmax,ymin,ymax,zmin,zmax)
      end
