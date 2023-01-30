**CHRDAT -- Convert seconds past 1 Jan 1985 to character
*+
      SUBROUTINE CHRDAT (SEC85, DATE)
      INTEGER       SEC85
      CHARACTER*(*) DATE
*
* Arguments:
* SEC85  (input) : Number of integral seconds past 1.0 Jan 1985
* DATE  (output) : Character containing date, hours, minutes, seconds in the
*                  form '910513 14:22:29'. Should be declared at least
*                  CHARACTER*15 in the calling (sub)program.
*-
* 13-May-1991: created - Remko Scharroo
*-----------------------------------------------------------------------
      INTEGER   TIME,MJD85,SEC,MIN,HOUR,DATUM,MDATE
      PARAMETER (MJD85=46066)

      TIME=SEC85
      SEC=MOD(TIME,60)
      TIME=TIME/60
      MIN=MOD(TIME,60)
      TIME=TIME/60
      HOUR=MOD(TIME,24)
      DATUM=MDATE(1,MJD85+TIME/24)
      WRITE (DATE,120) DATUM,HOUR,MIN,SEC

  120 FORMAT (I6.6,1X,I2.2,':',I2.2,':',I2.2)

      RETURN
      END
