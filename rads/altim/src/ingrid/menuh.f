      subroutine menuh(iopt,optio,quit)
*
* Prompt menu options
*
      include "ingrid.inc"
      character optio(MOPT)*41,frm*60,line*52,ans*41
      logical quit
      integer l1,l2,l3,iopt,nopt,i
*
 1000 format ('(',i2.2,'(''=''),1x,a',i2.2,',1x,',i2.2,'(''=''))')
 1010 format (/'InGrid-',i4,1x,a52,' used: ',i3,'-',i3,'%')
 1020 format (i1,'. ',a36,i2,'. ',a36)
 1030 format (i1,'. ',a36)
 1040 format (79('='))
 1050 format ('Enter choice -> ',$)
*
      quit=.false.
      l2=41
      ans=optio(1)
   10 if (ans(l2:l2).eq.' ') then
         l2=l2-1
         goto 10
      endif
      l1=(52-l2)/2
      l3=52-l1-l2
      write (frm,1000) l1-1,l2,l3-1
      write (line,frm) ans
      write (0,1010) MVRS,line,int(100./MA*jused+.5),
     .int(100./MA*(idnext-1)+.5)
      i=2
      do 15 nopt=-1,MOPT-3
   15    if (optio(nopt+3).eq.'-') goto 20
   20 do 25 i=0,nopt/2
         if (i+nopt/2+1.gt.nopt) then
            write (0,1030) i,optio(i+2)
         else
            write (0,1020) i,optio(i+2),i+nopt/2+1,optio(i+nopt/2+3)
         endif
   25 continue
      write (0,1040)
   30 write (0,1050)
      read (5,*,err=30,end=190) iopt
      if (iopt.lt.0 .or. iopt.gt.nopt) goto 30
      return
  190 quit=.true.
      end
