**LOWERCASE -- Change string to lower case
*+
      SUBROUTINE LOWERCASE (STRING)
      CHARACTER*(*) STRING
*
* This routine changes all upper case characters in a string to lower
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
          if (j.ge.65 .and. j.le.90) string(i:i)=char(j+32)
      enddo
      end
