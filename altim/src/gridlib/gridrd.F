#ifdef DOUBLE
#define GRIDRDn GRIDRD8
#define Real REAL*8
#define GRIDRDS GRIDRD8S
#define SUBNAME 'gridrd8'
#define BIG 1d35
#else
#define GRIDRDn GRIDRD4
#define Real REAL*4
#define GRIDRDS GRIDRD4S
#define SUBNAME 'gridrd4'
#define BIG 1e35
#endif
**GRIDRDn -- Read grid into a Real array
*+
      FUNCTION GRIDRDn (FILENM, NX, NY, Z,
     .                   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX)
      CHARACTER FILENM*(*)
      INTEGER*4 NX, NY, GRIDRDn
      Real    Z(*), XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
*
* This routine reads the contents of a grid file into a Real array.
* This array specified in the calling routine may be one- or two-dimensional.
* The values will be stored row by row, starting at the lower left corner,
* and pertain to the nodes of the grid. The coordinates of the value
* stored in Z(KX,KY) will thus be:
*
*   X = XMIN + (KX-1) * DX   with   DX = (XMAX-XMIN) / (NX-1)
*   Y = YMIN + (KY-1) * DY   with   DY = (YMAX-YMIN) / (NY-1)
*
* were NX and NY are the actual integer dimensions of the grid. The area
* covered by the grid is always rectangular and has corners (XMIN,YMIN) and
* (XMAX,YMAX).
* If the array is specified one-dimensional in the calling routine, you may
* specify an original X-dimension of 0 (NX=0) and NY equal to the maximum
* size of the single dimension. In that case, the array will be filled as
* were the original dimensions of the array equal to the actual dimensions
* of the grid.
* Points that do not have a defined value in the grid are given a value of
* 1.0D+35.
*
* Arguments:
*  FILENM   (input): Name of the file containing the grid.
*  NX       (input): X-dimension of array Z as defined in the calling
*                    routine. If NX=0, Z is defined one-dimensional.
*          (output): Actual X-dimension of the grid.
*  NY       (input): Y-dimension of array Z as defined in the calling routine.
*                    If NX=0, NY specifies the maximum single dimension.
*          (output): Actual Y-dimension of the grid, or zero if error occurred.
*  Z       (output): Array into which the gridded data must be stored.
*  XMIN    (output): X (longitude) value of lower left corner.
*  XMAX    (output): X (longitude) value of upper right corner.
*  YMIN    (output): Y (latitude) value of lower left corner.
*  YMAX    (output): Y (latitude) value of upper right corner.
*  ZMIN    (output): minimum Z (height) value in the grid.
*  ZMAX    (output): maximum Z (height) value in the grid.
*  GRIDRDn (output): Return code:
*                    =0 : no error, >0 : error
*-
* 11-Feb-1994 - created from former GRIDRD subroutine.
* 21-Jan-1997 - Improved opening of grid file. Better error trapping.
* 26-Nov-1997 - Include swap checking through compiler option -DSWAP.
* 20-Dec-2001 - Added extra argument to IOCONT because of change in FASTIO
* 17-Jun-2003 - Let FASTIO do stdin/stdout and remove sysdep.h
*-----------------------------------------------------------------------
      character line*80
      Real    zmid,fact
      integer*4 mx,my,maxi2,openf,readf,ios,fd,l
      parameter (maxi2=32767)

* Open grid file. Check if file exists.

      l=index(filenm,' ')-1
      if (l.le.0) l=len(filenm)
      fd=openf(filenm,'r')
      if (fd.lt.0) goto 1310

* Read header

      ios=readf(fd,80,line)
      if (ios.ne.80) goto 1320

* Interpret grid header

      if (line(1:5).ne.'@GRID') goto 1320
      read (line,600,err=1320) zmid,fact,mx,my,xmin,xmax,ymin,ymax
  600 format (5x,2e14.7,2i7,1x,4f8.3)
      zmin=zmid-maxi2*fact
      zmax=zmid+maxi2*fact
      if (xmin.eq.xmax.or.ymin.eq.ymax)
     |  write(0,1311)SUBNAME,'area not specified in ',filenm(:l)

* Check dimensions

      if (nx.eq.0) then
         if (mx*my.gt.ny) goto 1330
         nx=mx
      else
         if (mx.gt.nx .or. my.gt.ny) goto 1335
      endif

* Load entire data block.

      call GRIDRDS(fd,nx,ny,mx,my,zmid,fact,z,z)
      if (ny.eq.0) goto 1340
      nx=mx
      ny=my
      GRIDRDn=0
      goto 9999

* Error exits

 1310 write(0,1311)SUBNAME,'file not found: ',filenm(:l)
 1311 format (a,': ',a,a)
      GRIDRDn=2
      goto 1390

 1320 write(0,1311)SUBNAME,'illegal grid format in ',filenm(:l)
      GRIDRDn=2
      goto 1390

 1330 write(0,1331)SUBNAME,mx,my,ny,filenm(:l)
 1331 format (a,': grid dimension exceeds array size ',
     |        '(',i6,' x',i6,') > ',i9,' in ',a)
      GRIDRDn=3
      goto 1390

 1335 write(0,1336)SUBNAME,mx,my,nx,ny,filenm(:l)
 1336 format (a,': grid dimension exceeds array size ',
     |        '(',i6,' x',i6,') > (',i6,' x',i6,') in ',a)
      GRIDRDn=3
      goto 1390

 1340 write(0,1311)SUBNAME,'error loading grid ',filenm(:l)
      GRIDRDn=4
      goto 1390

 1390 ny=0

* Here normal termination

 9999 if (fd.ge.0) call closef(fd)
      end

* GRIDRDS is a subroutine that loads the INTEGER*2 image of the grid
* into the variable I and than converts it to Real values stored in Z.
* Note that I shares part of the memory of Z; therefore the transformation
* from I to Z is performed from back to front.

      subroutine GRIDRDS(fd,nx,ny,mx,my,zmid,fact,z,i)
      integer*4 fd,nx,ny,mx,my,kx,ky
      Real    z(nx,ny),zmid,fact,undefined
      integer*2 i(mx,my),j,mini2
      integer*4 readf,ios
      logical	ltlend
      parameter (undefined=BIG,mini2=-32768)

* Load all data at once

      ios=readf(fd,mx*my*2,i)
      if (ios.ne.mx*my*2) then
         ny=0
         return
      endif

* Swap integers if machine is Little Endian.

      if (ltlend()) call i2swap(mx*my,i)

* Process the data block. Convert integers to reals, back to front.
* Grid values -32768 are set to UNDEFINED.

      do ky=my,1,-1
         do kx=mx,1,-1
            j=i(kx,ky)
            if (j.eq.mini2) then
               z(kx,ky)=undefined
            else
               z(kx,ky)=zmid+fact*j
            endif
         enddo
      enddo
      end
