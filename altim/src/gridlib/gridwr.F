#ifdef DOUBLE
#define GRIDWRn GRIDWR8
#define Real REAL*8
#define SUBNAME 'gridwr8'
#define BIG 1d20
#else
#define GRIDWRn GRIDWR4
#define Real REAL*4
#define SUBNAME 'gridwr4'
#define BIG 1e20
#endif
**GRIDWRn -- Write Real array to a file
*+
      FUNCTION GRIDWRn (FILENM, NX, NY, Z, MX, XMIN, XMAX, YMIN, YMAX)
      INTEGER*4 NX, NY, MX, GRIDWRn
      Real    Z(MX,*), XMIN, XMAX, YMIN, YMAX
      CHARACTER FILENM*(*)
*
* This routine writes the contents of a Real array into a grid file.
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
* If the array is specified one-dimensional in the calling routine, you should
* specify MX=NX in the call to GRIDWRn.
* If values in the array are larger than 1D20, they will be labeled 'undefined'
* in the grid. Upon the next GRIDRD, a value of 1D35 will be returned for these
* points.
*
* Arguments:
*  FILENM   (input): Name of the file containing the grid.
*  NX       (input): Actual X-dimension of the grid. 
*  NY       (input): Actual Y-dimension of the grid.
*  Z        (input): Array in which the gridded data is stored.
*  MX       (input): X-dimension of the array as defined by the calling
*                    routine or program.
*  XMIN     (input): X (longitude) value of lower left corner.
*  XMAX     (input): X (longitude) value of upper right corner.
*  YMIN     (input): Y (latitude) value of lower left corner.
*  YMAX     (input): Y (latitude) value of upper right corner.
*  GRIDWRn (output): Return code:
*                    =0 : normal return, >0 : error occurred
*-
* 15-Dec-1993 - created for the previous gridwr subroutine.
* 26-Nov-1997 - Include swap checking though compiler option -DSWAP
* 17-Jun-2003 - Let FASTIO do stdin/stdout and remove sysdep.h
*-----------------------------------------------------------------------
      character line*80
      Real    fact,zmid,undefined,val,zmin,zmax
      integer*4 ndata,ipos
      parameter (ndata=1024)
      integer*4 openf,writef,fd,ios,mini2,kx,ky,l
      integer*2 jdata(ndata)
      logical	ltlend
      parameter (mini2=-32768,undefined=BIG)

      l=index(filenm,' ')
      fact=1
      if (mx.lt.1) goto 1310
      if (nx.gt.mx) goto 1320
      if (nx.lt.1 .or. ny.lt.1) goto 1330

* Determine extremes

      zmax=-undefined
      zmin=+undefined
      do ky=1,ny
         do kx=1,nx
            val=z(kx,ky)
            if (abs(val).lt.undefined) then
               zmin=min(val,zmin)
               zmax=max(val,zmax)
            endif
         enddo
      enddo

* Open file and write file header

      zmid=(zmin+zmax)/2
      if (zmax.ne.zmin) fact=(zmax-zmin)/65534
      write (line,600) zmid,fact,nx,ny,xmin,xmax,ymin,ymax
  600 format ('@GRID',2e14.7,2i7,1x,4f8.3)

      fd=openf(filenm,'w')
      if (fd.lt.0) goto 1340
      ios=writef(fd,80,line)
      if (ios.ne.80) goto 1350

* Write data

      ipos=0
      do ky=1,ny
         do kx=1,nx
            ipos=ipos+1
            if (ipos.gt.ndata) then
               if (ltlend()) call i2swap(ndata,jdata)
               ios=writef(fd,ndata*2,jdata)
               ipos=1
            endif
            val=z(kx,ky)
            if (abs(val).lt.undefined) then
               jdata(ipos)=nint((val-zmid)/fact)
            else
               jdata(ipos)=mini2
            endif
         enddo
      enddo
      if (ltlend()) call i2swap(ipos,jdata)
      ios=writef(fd,ipos*2,jdata)
      if (ios.ne.ipos*2) goto 1350
      GRIDWRn=0
      goto 9999

 1310 write(0,1311)SUBNAME,
     |'illegal x-dimension of matrix (mx<1)'
 1311 format (a,': ',a)
      GRIDWRn=1
      goto 9999

 1320 write(0,1311)SUBNAME,'x-dimension too large'
      GRIDWRn=2
      goto 9999

 1330 write(0,1311)SUBNAME,'illegal dimensions'
      GRIDWRn=3
      goto 9999

 1340 write(0,1341)
     |SUBNAME,'error opening output file ',filenm(:l)
 1341 format (a,': ',a,a)
      GRIDWRn=4
      goto 9999

 1350 write (0,1341)SUBNAME,'error writing grid buffer to ',
     |	filenm(:l)
      GRIDWRn=5
      goto 9999

* Here normal termination

 9999 call closef(fd)
      end
