**DAREAD -- Read sequentially from ADR file, blocked
*+
      SUBROUTINE DAREAD (IFILE, IREC, NRRECS, RDATA, ADATA)
      INTEGER IFILE, IREC, NRRECS
      REAL*8  RDATA(*)
      INTEGER*2 ADATA(*)
*
* Read record number IREC from ADR file number IFILE and store the
* content in RDATA. This reading is done from a block that is
* stored in memory. If a record beyond the block content is requested
* a new block is read from disk.
*
* Arguments
*   IFILE (input) : File number
*   IREC  (input) : Record number (starts with 1 for the first data record)
*   NRRECS (input) : Total number of data records in the file
*   RDATA (output) : Real data: time (sec), lat (rad), lon (rad),
*		     orbit (km), sea height (m)
*   ADATA (output) : Auxiliary data
*-
* 20-Oct-1994 - Bug: first measurement in block is wrong. Fixed
*               by adding extra line i=irec-blkrec(0)+1
*-----------------------------------------------------------------------
      integer*4 i
      include "maxcom.inc"

      i=irec-blkrec(0)+1      
      if (ifile.ne.blkfil(0) .or. i.gt.blkobs(0) .or. i.le.0)
     .   then
         blkrec(0)=irec
         blkobs(0)=min(nrrecs-blkrec(0)+1,maxobsblk)
         blkfil(0)=ifile
         call blload(1,0)
         i=irec-blkrec(0)+1      
      endif
      call blread(1,i,rdata,adata)
      end
