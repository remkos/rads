*************************************************************
      subroutine noclob(kb,error,quit)
*
* Check for unwritten changes
*
      implicit none
      include "ingrid.inc"
      character*4 answer
      logical error,quit
      integer kb
*
  550 format (a)
 1300 format ('ingrid: illegal buffer number:',i3)
 1310 format ('ingrid: Unwritten changes on buffer',i3,' will be lost'/
     .8x,'Continue? -> ',$)
*
      error=.true.
      quit=.false.
      if (kb.lt.1 .or. kb.gt.MBUF) then
         write (0,1300) kb
         kb=0
         return
      else if (used(kb) .and. fname(kb).eq.' ' .and. .not.proc(kb)) then
         write (0,1310) kb
         read (5,550,end=190,err=100) answer
         if (answer(1:1).eq.'n') return
      endif
      error=.false.
  100 return
  190 quit=.true.
      end
