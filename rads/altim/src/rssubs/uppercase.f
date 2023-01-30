**UPPERCASE -- Change string to upper case
*+
      SUBROUTINE UPPERCASE (STRING)
      CHARACTER*(*) STRING
*
* This routine changes all lower case characters in a string to upper
* case characters.
*
* Argument:
*  STRING (input/output): character string
*-
* 14-Nov-2001 : Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer i,j,n
      n=len(string)
      do i=1,n
          j=ichar(string(i:i))
          if (j.ge.97 .and. j.le.122) string(i:i)=char(j-32)
      enddo
      end
