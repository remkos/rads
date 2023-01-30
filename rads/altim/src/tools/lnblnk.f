**LNBLNK -- Determine the length of a string
*+
      FUNCTION LNBLNK (STRING)
      INTEGER*4 LNBLNK
      CHARACTER*(*) STRING
*
* This function returns the length of a character string, not counting
* trailing spaces. The function returns the full length of the character
* string as defined in the calling (sub)program if it ends in a character
* other than a space. If the string contains only spaces, a length of
* zero is returned.
*
* Arguments:
*  STRING  (input) : Character string of which the length is to be
*                    determined.
*  LNBLNK (output) : Length of the string (trailing spaces excluded)
*                    = 0, if the string contains only spaces
*                    = LEN(STRING), if the string ends in a non-space
*-
*  4-Mar-1992 : Created - Remko Scharroo
* 19-Mar-1992 : use I as loop variable, not LNBLNK
*-----------------------------------------------------------------------
      INTEGER*4 I
      DO I=LEN(STRING),1,-1
         IF (STRING(I:I).NE.' ') GOTO 90
      ENDDO
   90 LNBLNK=I
      END
