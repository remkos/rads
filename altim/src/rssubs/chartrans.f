**CHARTRANS -- Translate characters in a string
*+
      SUBROUTINE CHARTRANS (STRING, FROM, TO)
      CHARACTER*(*) STRING, FROM, TO

* This program translates each matching character of FROM in the string
* STRING with the corresponding character in TO.
*
* FROM and TO should be of the same length. If not, only the lesser
* number of characters are considered.
*
* Example:
*     STRING = 'original'
*     CALL CHARTRANS (STRING, 'oi', 'O!')
* results in:
*     STRING = 'Or!g!nal'
*-
* 13-Feb-2007 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer i,j,n

      n=min(len(from),len(to))

      do i=1,len(string)
         do j=1,n
	    if (string(i:i).eq.from(j:j)) then
	       string(i:i)=to(j:j)
	       goto 100
	    endif
	 enddo
100	 continue
      enddo
      end
