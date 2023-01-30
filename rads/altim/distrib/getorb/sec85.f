**SEC85 -- Convert MJD or YYMMDD or YYDDD to SEC85
*+
      FUNCTION SEC85 (I, DATE)
      REAL*8 SEC85, DATE
      INTEGER*4 I
*
* This function converts Modified Julian Dates (MJD) or Year-Month-Day
* (YYMMDD) or Year-Day (YYDDD) or Year-Month-Day-Hours-Minutes-Seconds
* (YYMMDDHHMMSS) to seconds from 1.0 Jan 1985 (SEC85). Fractions of days or
* seconds can also be included.
*
* Except for a 2-digit year indication (YY) it is now also possible to use
* a 4-digit indication (YYYY). The routine automatically recognises the
* form: whether YY or YYYY is specified.
*
* The first parameter (I) defines the input format: I=1 indicates MJD,
* I=2 means YYMMDD or YYYYMMDD, I=3 specifies YYDDD or YYYYDDD, and
* I=4 implies YYMMDDHHMMSS or YYYYMMDDHHMMSS, I=5 is a combination of I=2
* and I=4 (i.e, it takes YYMMDD, YYYYMMDD, YYMMDDHHMMSS and YYYYMMDDHHMMSS,
* as input).
*
* One can also choose NOT to specify the input format and let the function
* 'guess' which format it is. Obviously, there are limitations. it recognises
*         YYMMDD :  000101 - 040101  ( 2000/01/01 - 2004/01/01 )
*            MJD :   40587 - 53005   ( 1970/01/01 - 2004/01/01 )
*          YYDDD :   53005 - 99365   ( 1953/01/05 - 2000/01/01 )
*         YYMMDD :  500101 - 991232  ( 1950/01/01 - 2000/01/01 )
*        YYYYDDD :  1950d3 - 2050d3  ( 1950/01/01 - 2050/01/01 )
*       YYYYMMDD :  1950d4 - 2050d4  ( 1950/01/01 - 2050/01/01 )
*   YYMMDDHHMMSS :     1d8 - 1d12    ( 1950/01/01 - 2050/01/01 )
* YYYYMMDDHHMMSS : 1950d10 - 2050d10 ( 1950/01/01 - 2050/01/01 )
*
* Arguments:
*  I      (input): I=1, convert MJD to SEC85
*                  I=2, convert YYMMDD or YYYYMMDD to SEC85
*                  I=3, convert YYDDD or YYYYDDD to SEC85
*                  I=4, convert YYMMDDHHMMSS or YYYYMMDDHHMMSS to SEC85
*                  I=5, convert [YY]YYMMDD[HHMMSS] to SEC85
*                  I=0, convert any of the above to SEC85
*  DATE   (input): Input value.
*  SEC85 (output): Output value, seconds from 1.0 Jan 1985.
*-
*        1995 - Created. Inclusion of YYDDD and YYMMDDHHMMSS.
*  8-Jan-1996 - If/elseif's replaced by goto ().
*  5-Aug-1998 - YYYY notation added.
*  3-Nov-1999 - Add option I=5.
*-----------------------------------------------------------------------
      integer*4 yymmdd,mdate,mjd85,day,hh,mm
      real*8 ss
      parameter (mjd85=46066,day=86400)

      goto (1,10,20,30,40,50) i+1
      SEC85=date
      return

1     if (date.lt.0) then
	 SEC85=date
         return
      else if (date.lt.40587) then
         goto 20
      else if (date.lt.53005) then
         goto 10
      else if (date.lt.500101) then
         goto 30
      else if (date.lt.1950d3) then
         goto 20
      else if (date.lt.1950d4) then
	 goto 30
      else if (date.le.1d8) then
	 goto 20
      else if (date.lt.2050d10) then
	 goto 40
      else
	 SEC85=date
         return
      endif

10    continue
* MJD -> SEC85
	 sec85=(date-mjd85)*day
	 return
20    continue
* YYMMDD or YYYYMMDD -> SEC85
	 yymmdd=int(date)
	 sec85=(mdate(2,yymmdd)-mjd85+date-yymmdd)*day
	 return
30    continue
* YYDDD or YYYYDDD -> SEC85
	 yymmdd=int(date/1d3)*10000+0101
	 sec85=(mdate(2,yymmdd)-mjd85+mod(date,1000d0)-1)*day
	 return
40    continue
* YYMMDDHHMMSS or YYYYMMDDHHMMSS -> SEC85
	 yymmdd=int(date/1d6)
	 ss=mod(date,1d6)
	 hh=int(ss/1d4)
	 ss=ss-hh*1d4
	 mm=int(ss/1d2)
	 ss=ss-mm*1d2
	 sec85=(mdate(2,yymmdd)-mjd85)*day+hh*3600+mm*60+ss
	 return
50    continue
* YYMMDD or YYYYMMDD or YYMMDDHHMMSS or YYYYMMDDHHMMSS -> SEC85
         if (date.lt.1d8) then
            goto 20
         else
            goto 40
         endif
      end
