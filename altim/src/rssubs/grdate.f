**GRDATE -- Get current date
*+
      SUBROUTINE GRDATE (DATE)
      CHARACTER*(*) DATE

* This routine returns the current date in format '20-Jan-1998 18:22:13'
* (this is a total of 20 characters, the rest of DATE will be filled with
* spaces).
*
* Output argument:
*  DATE : The current date in the format 20-Jan-1998 18:22:13
*-
*  4-Nov-1998 -- Created by Remko Scharroo
*-----------------------------------------------------------------------
      character temp*26
      call fdate(temp)
      date=' '
      date( 1: 2)=temp( 9:10)
      date( 3: 3)='-'
      date( 4: 6)=temp( 5: 7)
      date( 7: 7)='-'
      date( 8:11)=temp(21:24)
      date(12:20)=temp(11:19)
      end
