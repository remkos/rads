      subroutine newbuf(text,kb,nxp,nyp,error,quit)
*
* Ask for destination buffer, check if this buffer is new or old and
* check if there is enough space within the buffer to contain the
* predefined dimensions nxp and nyp. If not shift the starting position
* of the buffer.
*
      include "ingrid.inc"
      character*(*) text
      integer kb,nxp,nyp,nspace,idmin,i,j,lb
      logical error,quit,suffic
*
  550 format (a)
  551 format (a,$)
*
* Ask for target buffer
*
      error=.true.
      quit=.false.
      write (0,551) text
      write (0,551) ' -> '
      read (5,*,end=190,err=100) kb
      call noneg(kb)
      if (kb.lt.1 .or. kb.gt.MBUF) then
         write (0,550) 'ingrid: buffer number out of range'
         return
      else if (blkd(kb)) then
         write (0,550)
     .'ingrid: cannot use buffer more than once during this action'
         return
      endif
*
* If buffer is used, determine if buffer has sufficient space
*
      nspace=nxp*nyp
      suffic=.false.
      if (used(kb)) then
         idmin=MA+1
         j=id(kb)+nx(kb)*ny(kb)
         do 20 lb=1,MBUF
            i=id(lb)
   20       if (used(lb) .and. i.ge.j .and. i.lt.idmin) idmin=i
         suffic=(id(kb)+nspace.le.idmin)
      endif
*
      if (.not.suffic) then
*
* If (1) buffer is not used or (2) buffer is too small,
* create buffer elsewhere
*
         if (idnext+nspace.gt.MA+1) then
            write (0,550) 'ingrid: No more space; '//
     .      'wipe some buffers first'
         else
*
* Check for unwritten changes
*
            call noclob(kb,error,quit)
            if (quit) goto 190
            if (error) return
*
* Start buffer at end of array a
*
            id(kb)=idnext
            idnext=idnext+nspace
            nx(kb)=nxp
            ny(kb)=nyp
            error=.false.
         endif
      else
*
* Grid will fit in used buffer; check for unwritten changes
*
         call noclob(kb,error,quit)
         if (quit) goto 190
         if (error) return
         nx(kb)=nxp
         ny(kb)=nyp
         error=.false.
      endif
      proc(kb)=.false.
      fname(kb)=' '
  100 return
  190 quit=.true.
      end
