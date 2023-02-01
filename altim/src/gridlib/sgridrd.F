#ifdef DOUBLE
#define SGRIDRD SGRIDRD8
#define Real REAL*8
#define SGRIDRS SGRIDRD8S
#define SUBNAME 'sgridrd8'
#define BIG 1d35
#else
#define SGRIDRD SGRIDRD4
#define Real REAL*4
#define SGRIDRS SGRIDRD4S
#define SUBNAME 'sgridrd4'
#define BIG 1e35
#endif
**SGRIDRD -- Read grid into a Real array
*+
      FUNCTION SGRIDRD (FILENM, NX, NY, Z,
     .                   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX)
      CHARACTER FILENM*(*)
      INTEGER*4 NX, NY, SGRIDRD
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
* Subgrids.
* If the grid boudaries are not given a priori (i.e. when XMIN=XMAX or
* YMIN=YMAX), then the entire grid as it is stored in returned, along with
* its actual boundaries. If the boundaries are predefined, a subsection of
* the grid is selected. Still, the actual boundaries are returned, which
* may differ from the a priori because the grid mesh may not comply with
* the selected boundaries.
*
* Sampling.
* If the array size declared in the calling routine (defined by NX and NY)
* is too small to contain the grid (or the requested subsection), the
* grid is sampled. Only one out of 1, 2, 3, ... grid points in each
* direction will be returned.
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
*  XMIN    (in/output): X (longitude) value of lower left corner.
*  XMAX    (in/output): X (longitude) value of upper right corner.
*  YMIN    (in/output): Y (latitude) value of lower left corner.
*  YMAX    (in/output): Y (latitude) value of upper right corner.
*  ZMIN    (output): minimum Z (height) value in the grid.
*  ZMAX    (output): maximum Z (height) value in the grid.
*  SGRIDRD (output): Return code:
*                    =0 : no error, >0 : error
*-
*    Jan-1995 : Created by Remko Scharroo
* 22-Feb-1995 : Change selection boundaries. Just give a little bit more
*               than requested. This is to have one cell extra beyond the
*               boundary so interpolation is possible along the boundary.
* 24-Feb-1995 : Manual update. Adjusted vertical boundary selection.
* 28-Oct-1996 : Include wrapping at BOTH sides.
* 29-Sep-2001 : Improved sizing of subgrids (esp. when size is even)
*-----------------------------------------------------------------------
#include <sysdep.h>
      character line*80
      Real    zmid,fact,x0,x1,y0,y1
      real*8  dx,dy
      integer*4 gx,gy,mx,my,mx0,mx1,my0,my1,maxi2,xscale,yscale,
     |          openf,readf,ios,fd,l
      logical   exist,info

      parameter (maxi2=32767)

* Open grid file. First check if file exists.

      l=index(filenm,' ')-1
      if (l.le.0) l=len(filenm)
      inquire (exist=exist,file=filenm)
      if (.not.exist) goto 1310

* Read header

      fd=openf(filenm,'r')
      ios=readf(fd,80,line)
      if (ios.ne.80) goto 1320

* Interpret grid header

      if (line(1:5).ne.'@GRID') goto 1320
      read (line,600,err=1320) zmid,fact,mx,my,x0,x1,y0,y1
  600 format (5x,2e14.7,2i7,1x,4f8.3)
      zmin=zmid-maxi2*fact
      zmax=zmid+maxi2*fact

* Use the standard GRIDRD routine:
* - If area not specified on call
* - If area not specified in grid

      info=(ny.lt.0)
      ny=abs(ny)

      if (x0.eq.x1 .or. y0.eq.y1) then
         x0=xmin
         x1=xmax
         y0=ymin
         y1=ymax
      else if (xmin.eq.xmax .or. ymin.eq.ymax) then
         xmin=x0
         xmax=x1
         ymin=y0
         ymax=y1
      endif

* Determine what block to load. Just a bit larger if you can, please.

      dx=(x1-x0)/(mx-1)
      dy=(y1-y0)/(my-1)
      mx0=nint((xmin-x0)/dx+0.501d0)
      mx1=nint((xmax-x0)/dx+1.499d0)
      if (mx1.gt.mx .and. x1-x0.lt.359) mx1=mx
      if (mx1.lt. 1 .and. x1-x0.lt.359) mx0=1
      my0=nint((ymin-y0)/dy+0.501d0)
      my1=nint((ymax-y0)/dy+1.499d0)
      if (my0.lt. 1) my0=1
      if (my1.gt.my) my1=my

* Set the apriori size of the subgrid. Start with scale=1.

      gx=mx1-mx0+1
      gy=my1-my0+1
      xscale=1
      yscale=1

* Check dimensions.

      if (nx.ne.0) then

* If the array in the calling routine is two-dimensional reduce both
* grid dimension to below the specified dimensions.

  110    gx=(mx1-mx0)/xscale+1
         if (gx.gt.nx) then
            xscale=xscale+1
            goto 110
         endif
  112    gy=(my1-my0)/yscale+1
         if (gy.gt.ny) then
            yscale=yscale+1
            goto 112
         endif

      else

* If the array is one-dimensional, decrease first the one that has been
* scaled less, or (if scaled equally) decrease the larger dimension.

  120    gx=(mx1-mx0)/xscale+1
         gy=(my1-my0)/yscale+1
         if (gx*gy.gt.ny) then
            if (xscale.lt.yscale) then
               xscale=xscale+1
            else if (xscale.gt.yscale) then
               yscale=yscale+1
            else if (gx.gt.gy) then
               xscale=xscale+1
            else
               yscale=yscale+1
            endif
            goto 120
         endif

      endif

* Scales are now determined. Reset the grid dimensions.
* If they differ from the original, give a printout.

      mx1=mx0+(gx-1)*xscale
      my1=my0+(gy-1)*yscale
      xmin=nint((x0+(mx0-1)*dx)*1d6)/1d6
      xmax=nint((x0+(mx1-1)*dx)*1d6)/1d6
      ymin=nint((y0+(my0-1)*dy)*1d6)/1d6
      ymax=nint((y0+(my1-1)*dy)*1d6)/1d6
      if (nx.eq.0) nx=gx

      if ((gx.ne.mx .or. gy.ne.my) .and. .not.info) then
         write (0,1311)
     |		SUBNAME,'Selecting subgrid of file: ',filenm(:l)
         write (0,1301) '... longitude',
     |		x0,x1,mx,xmin,xmax,mx0,mx1,xscale,gx
         write (0,1301) '... latitude ',
     |		y0,y1,my,ymin,ymax,my0,my1,yscale,gy
 1301    format (a,':',2f8.2,i6,' ->',2f8.2,2i6,i4,i6)
      endif

* Load selected data block. We have:
* X0, X1 : grid longitude boundaries
* Y0, Y1 : grid latitude boundaries
* XMIN, XMAX : selected longitude boundaries
* YMIN, YMAX : selected latitude boundaries
* NX, NY : the dimensions of the array as defined in the calling routine
* MX, MY : the dimensions of the grid as stored
* GX, GY : the dimensions of the selected grid
* MX0, MX1, XSCALE : the first and last grid column to be read with skips
* MY0, MY1, YSCALE : the first and last grid line to be read with skips

      if (info) then
      else if (mx0.lt.1 .and. mx1.gt.mx) then

* If there is data to be wrapped from the rightside of the grid to the left,
* AND from the leftside to the right.
* 1. Read the middle and store it there
* 2. Grab the lefthand side and store it right
* 3. Grab the righthand side and store it left

         zmin=1e35
         zmax=-1e35
         call SGRIDRS(fd,nx,ny,mx,1,mx,xscale,my0,my1,yscale,
     |      zmid,fact,z((1-mx0)/xscale+1),zmin,zmax)
                if (ny.eq.0) goto 1340
         mx0=max(nint((xmin-360-x0)/dx+1),1)
         mx1=nint((xmax-360-x0)/dx+1)
         call SGRIDRS(fd,nx,ny,mx,mx0,mx1,xscale,my0,my1,yscale,
     |      zmid,fact,z(gx-(mx1-mx0)/xscale),zmin,zmax)
         mx0=nint((xmin+360-x0)/dx+1)
         mx1=min(nint((xmax+360-x0)/dx+1),mx)
         call SGRIDRS(fd,nx,ny,mx,mx0,mx1,xscale,my0,my1,yscale,
     |      zmid,fact,z(1),zmin,zmax)

      else if (mx0.lt.1) then

* If there is data to be wrapped from the rightside of the grid to the left,
* first read the rightside and store it left, then the leftside and store
* it right

         zmin=1e35
         zmax=-1e35
         call SGRIDRS(fd,nx,ny,mx,1,mx1,xscale,my0,my1,yscale,
     |      zmid,fact,z(gx-(mx1-1)/xscale),zmin,zmax)
         if (ny.eq.0) goto 1340
         mx0=nint((xmin+360-x0)/dx+1)
         mx1=min(nint((xmax+360-x0)/dx+1),mx)
         call SGRIDRS(fd,nx,ny,mx,mx0,mx1,xscale,my0,my1,yscale,
     |      zmid,fact,z(1),zmin,zmax)

      else if (mx1.gt.mx) then

* If there is data to be wrapped from the leftside of the grid to the right,
* first read the rightside and store it left, then the leftside and store
* it right

         zmin=1e35
         zmax=-1e35
         call SGRIDRS(fd,nx,ny,mx,mx0,mx,xscale,my0,my1,yscale,
     |      zmid,fact,z(1),zmin,zmax)
         if (ny.eq.0) goto 1340
         mx0=max(nint((xmin-360-x0)/dx+1),1)
         mx1=nint((xmax-360-x0)/dx+1)
         call SGRIDRS(fd,nx,ny,mx,mx0,mx1,xscale,my0,my1,yscale,
     |      zmid,fact,z(gx-(mx1-mx0)/xscale),zmin,zmax)


      else

* Just read the selection and store it

         zmin=1e35
         zmax=-1e35
         call SGRIDRS(fd,nx,ny,mx,mx0,mx1,xscale,my0,my1,yscale,
     |      zmid,fact,z(1),zmin,zmax)
      endif

      nx=gx
      ny=gy
      SGRIDRD=0
      goto 9999

* Error exits

 1310 write (0,1311) SUBNAME,'file not found: ',filenm(:l)
 1311 format (a,': ',a,a)
      SGRIDRD=2
      goto 1390

 1320 write(0,1311) SUBNAME,'illegal grid format in ',filenm(:l)
      SGRIDRD=2
      goto 1390

 1340 write(0,1311) SUBNAME,'error loading grid ',filenm(:l)
      SGRIDRD=4
      goto 1390

 1390 ny=0

* Here normal termination

 9999 call closef(fd)
      end

* SGRIDRS is a subroutine that loads the INTEGER*2 image of the grid
* into the variable I and than converts it to Real values stored in Z.
* Note that I shares part of the memory of Z; therefore the transformation
* from I to Z is performed from back to front.

      subroutine SGRIDRS(fd,nx,ny,mx,mx0,mx1,dx,my0,my1,dy,
     |		zmid,fact,z,zmin,zmax)
      integer*4 fd,nx,ny,mx,gx,gy,mx0,mx1,dx,my0,my1,dy
      integer*4 readf,seekf,ios,kx,ky,ix,iy,bufsize
      Real    z(nx,ny),zmid,fact,undefined,zmin,zmax
      parameter (undefined=BIG,bufsize=25000)
      integer*2 buf(bufsize),j,mini2
      parameter (mini2=-32768)

      gx=mx1-mx0+1
      gy=my1-my0+1

* Check if buffer size is sufficient

      if (gx*2.gt.bufsize) then
         write (0,
     |		"('SGRIDRD: Buffer size insufficient; need ',i6)") gx*2
         call fin("Increase buffer size or reduce grid")
      endif

* Read all lines one by one

      ky=0
      do iy=my0,my1,dy
         ky=ky+1
         ios=seekf(fd,80+((iy-1)*mx+(mx0-1))*2,0)
         ios=readf(fd,gx*2,buf)

#ifdef SWAP
* Swap integers if machine is Little Endian.

         call i2swap(gx,buf)
#endif

* Process the data block. Convert integers to reals.
* Grid values -32768 are set to UNDEFINED.

         kx=0
         do ix=1,gx,dx
            kx=kx+1
            j=buf(ix)
            if (j.eq.mini2) then
               z(kx,ky)=undefined
            else
               z(kx,ky)=zmid+fact*j
               zmin=min(zmin,z(kx,ky))
               zmax=max(zmax,z(kx,ky))
            endif
         enddo
      enddo
      end
