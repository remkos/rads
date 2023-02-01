**CHECKENV -- Get environment variable or leave default
*+
      SUBROUTINE CHECKENV (ENV, STRING, L)
      CHARACTER*(*) ENV, STRING
      INTEGER*4     L

* This routine returns the contents of environment variable ENV in the
* variable STRING.
* If the environment variable ENV is not set, STRING is unchanged.
* Upon return the L will be the length of the string STRING, without
* trailing spaces.
*
* Arguments:
*   ENV     (input) : Name of the environment variable.
*   STRING  (input) : Default value for STRING.
*          (output) : Contents of the environment variable or default.
*   L      (output) : Length of STRING
*-
* 13-Feb-1998 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer lnblnk
      character*160 temp

      call getenv (ENV, temp)
      if (temp.ne.' ') then
         STRING=temp
      else if (STRING.eq.' ') then
	 l=lnblnk(ENV)
         write (0,600) ENV(:l)
      endif
      L=lnblnk(STRING)

600   format ('WARNING: No environment variable of default for ',a)
      end
