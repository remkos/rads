**LISTARGS -- List all input arguments
*+
      SUBROUTINE LISTARGS (UNIT, LENGTH)
      INTEGER UNIT, LENGTH
*
* This subroutine lists all input arguments from the command line,
* including the name of the command. Output will be written
* to unit number UNIT. Lines will be cut before the maximum length
* of the output lines, specified by LENGTH is encountered.
* Specify LENGTH=0 for infinitely long lines.
*
* Arguments:
*  UNIT   (input): Unit number for output
*  LENGTH (input): Maximum length of output lines (0 for infinite lines)
*-
*  3-Apr-1996 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4    i,l,iarg,iargc
      character*80 argum
      
      l=0
      do iarg=0,iargc()
         call getarg(iarg,argum)
         i=index(argum,' ')
         l=l+i
         if (l.ge.length .and. length.ne.0) then
            write (unit,550)
            l=i
         endif
         write (unit,551) argum(:i)
      enddo
      write (unit,550)

550   format(a)
551   format(a,$)

      end
