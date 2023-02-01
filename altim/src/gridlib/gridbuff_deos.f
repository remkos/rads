**GRIDBUFF_DEOS -- Load a grid into memory (DEOS format)
*+
      FUNCTION GRIDBUFF_DEOS (FILENM, POINTER)
      INTEGER*4 GRIDBUFF_DEOS, POINTER
      CHARACTER FILENM*(*)
*-
* 16-Aug-2006 - Introduced new OPENF
* 26-Jul-2006 - Avoiding use of %val
* 21-Aug-2005 - Created by Remko Scharroo from GRIDBUFF
*-----------------------------------------------------------------------
      character line*80
      integer*4 fd,ios,nbytes,openf,readf,mallocf,pntr2,nb,memloc
      logical	ltlend,lendian/.false./
      include "gridbuff.inc"

* Open grid file and read header

      fd=openf(filenm,'r')
      if (fd.lt.0) goto 1310
      ios=readf(fd,80,line)
      if (ios.ne.80) goto 1320

* Interpret grid header

      if (line(:3).ne.'@GR') goto 1350
      if (line(4:5).eq.'ID') line(4:5)='2B'
      read (line(4:4),*) nb
      if (nb.eq.1 .or. nb.eq.4) then
         ntype=nb
      else if (nb.eq.2) then
         ntype=3
      else
         goto 1320
      endif
      if (line(5:5).eq.'L') then
         lendian=.true.
      else if (line(5:5).eq.'B') then
         lendian=.false.
      else
         goto 1320
      endif
      read (line,600,err=1320) z0,dz,nx,ny,xmin,xmax,ymin,ymax
600   format (5x,2d14.7,2i7,1x,4f8.3)

* Determine dimensions and allocate memory

      nbytes=nx*ny*nb		! Number of bytes of grid data
      nbuf=(nbytes+mhead+7)/8*8	! Round up to multiple of 8 bytes
      if (mallocf(nbuf,pointer).gt.0) goto 1330
      pntr2=pointer+mhead-memloc(tmp_b)

* Store grid information in the data buffer header

      dx=(xmax-xmin)/(nx-1)
      dy=(ymax-ymin)/(ny-1)
      if (nb.eq.1) then
         zmin=z0-127*dz
         zmax=z0+127*dz
	 znan=-128
      else if (nb.eq.2) then
         zmin=z0-32767*dz
         zmax=z0+32767*dz
	 znan=-32768
      else
         zmin=z0-2147483647*dz
         zmax=z0+2147483647*dz
         znan=-2147483647-1
      endif
      call memput(pointer,mhead,head)

* Load entire data block.

      ios=readf(fd,nbytes,tmp_b(pntr2))
      if (ios.ne.nbytes) goto 1340

* Swap integers if machine`s endianness is not the same as the grid`s

      if (ltlend().eqv.lendian) then
      else if (nb.eq.4) then
         call i4swap(nx*ny,tmp_b(pntr2))
      else
         call i2swap(nx*ny,tmp_b(pntr2))
      endif
      GRIDBUFF_DEOS=0
      goto 9999

* Error exits

1310  GRIDBUFF_DEOS=1
      goto 9998
1320  GRIDBUFF_DEOS=2
      goto 9998
1330  GRIDBUFF_DEOS=3
      goto 9998
1340  GRIDBUFF_DEOS=4
      goto 9998
1350  GRIDBUFF_DEOS=5

* Here normal termination

9998  POINTER=0
9999  if (fd.ge.0) call closef(fd)
      end
