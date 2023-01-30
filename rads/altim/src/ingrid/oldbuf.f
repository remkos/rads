*************************************************************
      subroutine oldbuf(text,kb,negat,error,quit)
*
* Prompt for buffer with existing grid
*
      include "ingrid.inc"
      character*(*) text
      integer kb,kbtest
      logical error,quit,negat
*
  550 format (a)
  551 format (a,$)
*
      error=.true.
      quit=.false.
      write (0,551) text
      write (0,551) ' -> '
      read (5,*,end=190,err=100) kb
      if (.not.negat) call noneg(kb)
      kbtest=iabs(kb)
      if (kbtest.lt.1 .or. kbtest.gt.MBUF) then
        write (0,550) 'ingrid: illegal buffer number'
        return
      else if (.not.used(kbtest)) then
        write (0,550) 'ingrid: buffer not in use'
        return
      endif
      error=.false.
  100 return
  190 quit=.true.
      end
