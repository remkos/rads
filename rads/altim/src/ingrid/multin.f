      subroutine multin(text,mbuf,n,idata,error,quit)
**********************************************************************
* Read a number (max=mbuf) of integer values from standard input;    *
* return their values (idata) and amount (n)                         *
*                                                                    *
* Remko Scharroo, Delft, August 28, 1990.                            *
**********************************************************************
      integer idata(mbuf)
      logical error,quit,dash
      character text*(*),line*160,a*1,substr*160
      integer n,i0,last,i,j,mbuf
*
  550 format (a)
  551 format (a,$)
  560 format (i80)
*
      n=0
      i0=1
      last=0
      dash=.false.
      error=.true.
      quit=.false.
      write (0,551) text
      write (0,551) ' -> '
      read (5,550,end=190,err=100) line
      do 50 i=1,160
         a=line(i:i)
         if (a.eq.' ' .or. a.eq.',' .or. a.eq.'-') then
            if (i.gt.i0) then
               substr=line(i0:i-1)
               read (substr,560,err=40) idata(n+1)
               if (dash) then
                  do 30 j=last+1,idata(n+1)
                     n=n+1
                     idata(n)=j
                     if (n.ge.mbuf) goto 90
   30             continue
                  dash=.false.
               else
                  n=n+1
                  if (n.ge.mbuf) goto 90
               endif
               last=idata(n)
            endif
            if (a.eq.'-') dash=.true.
   40       i0=i+1
         endif
   50 continue
      if (dash) then
         do 60 i=1,mbuf-n
   60       idata(n+i)=last+i
         n=mbuf
      endif
   90 if (n.ge.1) error=.false.
  100 return
*
  190 quit=.true.
      end
