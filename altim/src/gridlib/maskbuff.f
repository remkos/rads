**MASKBUFF -- Load a mask file into memory
*+
      FUNCTION MASKBUFF (FILENM, POINTER)
      INTEGER*4 MASKBUFF, POINTER
      CHARACTER*(*) FILENM
*
* This routine allocated memory and loads the contents of a mask
* file into the allocated memory.
*
* The mask file has to be coded in the PBM format. In the
* header of the PBM file the mask area can be coded using a line
* # AREA <lon0> <lon1> <lat0> <lat1>
* If this line is absent the default will be assumed:
* # AREA 0 360 -90 90
*
* The routine MASKBINT can be used to interogate the mask buffer
* to determine the flag bit at a given latitude and longitude.
*
* Other information on the mask, like mask dimension, can be obtained
* with the MASKBINF routine.
*
* When the allocation of memory or the loading of the grid was
* unsuccessful, this will be reflected in the returned pointer value.
*
* Input argument:
*  FILENM   : Name of the BPM file containing the mask.
*
* Output argument:
*  POINTER  : Pointer to the mask structure
*
* Return value:
*  MASKBUFF : 0 = No error
*             1 = Could not find or open mask
*             2 = Ilegal mask format
*             3 = Memory allocation was unsuccessful
*             4 = Error loading mask
*-
* 26-Jul-2006 - Avoiding use of %val
* 22-Jul-2005 - Use memory allocation
* 28-Jan-2002 - Created from GRIDBUFF
*-----------------------------------------------------------------------
      character line*256
      integer*4 fd,ios,l,lnblnk,nbytes,readf,openf,
     |		mallocf,pntr2,memloc
      logical	ltlend
      include "gridbuff.inc"

* Open mask file. Check if file exists.

      l=lnblnk(filenm)
      fd=openf(filenm,'r')
      if (fd.lt.0) goto 1310

* Read header

      ios=readf(fd,0,line)-1
      if (line(:ios).ne.'P4') goto 1320

* Interpret mask header

      xmin=0d0
      xmax=360d0
      ymin=-90d0
      ymax=90d0
10    ios=readf(fd,0,line)-1
      if (line(:6).eq.'# AREA') then
         read (line(7:ios),*) xmin,xmax,ymin,ymax
         goto 10
      else if (line(:1).eq.'#') then
         goto 10
      else
         read (line(:ios),*) nx,ny
      endif

* Determine dimensions and allocate memory

      nbytes=(nx*ny+7)/8	! Total number of bytes of data in mask file
      nbuf=(nbytes+mhead+7)/8*8	! Round up to multiple of 8 bytes
      if (mallocf(nbuf,pointer).gt.0) goto 1330
      pntr2=pointer+mhead-memloc(tmp_b)

* Store mask information in header of the data buffer

      ntype=0
      z0=0
      dz=0
      zmin=0
      zmax=1
      dx=(xmax-xmin)/nx
      dy=(ymax-ymin)/ny
      call memput(pointer,mhead,head)

* Load entire data block.

      ios=readf(fd,nbytes,tmp_b(pntr2))
      if (ios.ne.nbytes) goto 1340

* Swap integers if machine is Little Endian.

      if (ltlend()) call i4swap((nbytes+3)/4,tmp_b(pntr2))
      MASKBUFF=0
      goto 9999

* Error exits

1310  write(0,1311) 'MASKBUFF: file not found: ',filenm(:l)
1311  format (a,a)
      MASKBUFF=1
      goto 9998

1320  write(0,1311) 'MASKBUFF: illegal mask format in ',filenm(:l)
      MASKBUFF=2
      goto 9998

1330  write(0,1311) 'MASKBUFF: error allocating memory for ',filenm(:l)
      MASKBUFF=3
      goto 9998

1340  write(0,1311) 'MASKBUFF: error loading mask ',filenm(:l)
      MASKBUFF=4

* Here normal termination

9998  POINTER=0
9999  if (fd.ge.0) call closef(fd)
      end
