**YMDHMS -- Convert seconds past 1 Jan 1985 to YYMMDDHHMMSS.SS
*+
      FUNCTION YMDHMS (SEC85)
      REAL*8 YMDHMS, SEC85
*
* Arguments:
* SEC     (input) : Number of seconds past 1.0 Jan 1985
* YMDHMS (output) : Real number formatted as YYMMDDHHMMSS.SS
*-
*  6-Aug-1997 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      INTEGER   TIME,MJD85,SEC,MIN,HOUR,DATUM,MDATE
      REAL*8    FRAC
      PARAMETER (MJD85=46066)

      TIME=INT(SEC85)
      FRAC=SEC85-TIME
      SEC=MOD(TIME,60)
      TIME=TIME/60
      MIN=MOD(TIME,60)
      TIME=TIME/60
      HOUR=MOD(TIME,24)
      DATUM=MDATE(1,MJD85+TIME/24)
      YMDHMS=DATUM*1D6+HOUR*1D4+MIN*1D2+SEC*1D0+FRAC

      RETURN
      END
