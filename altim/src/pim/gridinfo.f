**GRIDINFO -- Read grid into a REAL*4 array
*+
      FUNCTION GRIDINFO (FILENM, NX, NY,
     .                   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX)
      CHARACTER FILENM*(*)
      INTEGER*4 NX, NY, GRIDINFO
      REAL*4    XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
*
* This routine reads the information of a grid file.
*
* Arguments:
*  FILENM   (input): Name of the file containing the grid.
*  NX      (output): Actual X-dimension of the grid.
*  NY      (output): Actual Y-dimension of the grid, or zero if error occurred.
*  Z       (output): Array into which the gridded data must be stored.
*  XMIN    (output): X (longitude) value of lower left corner.
*  XMAX    (output): X (longitude) value of upper right corner.
*  YMIN    (output): Y (latitude) value of lower left corner.
*  YMAX    (output): Y (latitude) value of upper right corner.
*  ZMIN    (output): minimum Z (height) value in the grid.
*  ZMAX    (output): maximum Z (height) value in the grid.
*  GRIDINFO(output): Return code:
*                    =0 : no error, >0 : error
*-----------------------------------------------------------------------
      integer sgridrd4,l1,l2
      real dum

      l1=index(filenm,'{')+2
      if (l1.lt.6) l1=1
      l2=index(filenm(l1:),' ')+l1-2
      if (filenm(l1:l1).eq.'(' .and. filenm(l2:l2).eq.')') then
	 nx=2
	 ny=2
	 gridinfo=0
      else
         ny=-ny
         gridinfo=sgridrd4(filenm(l1:l2),nx,ny,
     &		dum,xmin,xmax,ymin,ymax,zmin,zmax)
      endif
      end
