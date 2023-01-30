**MDATE -- Convert MJD to (YY)YYMMDD or v.v.
*+
      FUNCTION MDATE (I, J)
      implicit none
      INTEGER  MDATE, I, J
*
* This function converts Modified Julian Dates (MJD) to Year-Month-Day
* in the format YYMMDD or YYYYMMDD, or vice versa.
*
* Dates in the 21st century are given in 6-digits, similar to dates in
* the 20th century. So 13 October 2001 is written as 011013.
* Use I=1 to convert MJD to YYMMDD and I=2 to convert YYMMDD to MJD.
* The algorithm is valid for the years 1955 till 2054 only.
*
* Optionally, one can enter the YYYYMMDD including the century indication,
* or have it returned by MDATE. Thus 13 October 2001 is written as 20011013.
* Use I=2 to convert YYYYMMDD to MJD and I=3 to convert MJD to YYYYMMDD.
* With these options, the algorithm is valid for 1 March 1900 until
* 28 February 2100.
*
* When YYYYMMDD or YYMMDD is input to the function, a month number of
* 0 or 13 will also work. Day numbers 0 and 32 are also excepted.
* For example: 000000 will be recognized as 30 November 1999.
*
* Arguments:
*  I      (input): I=1, convert        MJD -> YYMMDD
*                  I=2, convert (YY)YYMMDD -> MJD
*                  I=3, convert        MJD -> YYYYMMDD
*  J      (input): Either MJD, YYMMDD, or YYYYMMDD (depending on I).
*  MDATE (output): Either YYMMDD, MJD, or YYYYMMDD (depending on I).
*-
* 15-Sep-1988 - Created: Remko Scharroo, DUT/SOM (c), DUT/KOSG
* 14-Jan-1992 - Revised version
*  6-Aug-1998 - Limited range to 1955-2054 instead of 1960-2059. Added
*               YYYYMMDD input and output
* 21-May-2004 - Added posibility to use month=0.
*-----------------------------------------------------------------------
      integer   t,y,m,d

      if (i.eq.2) then
          d=mod(j,100)
          t=(j-d)/100
          m=mod(t,100)
          y=(t-m)/100
	  call ymd2mjd(y,m,d,mdate)
      else
	  call mjd2ymd(j,y,m,d)
	  if (i.eq.1) y=mod(y,100)
          mdate=y*10000+m*100+d
      endif
      end

**YMD2MJD -- Convert YY, MM and DD to MJD
*+
      SUBROUTINE YMD2MJD (YY, MM, DD, MJD)
      implicit none
      INTEGER YY, MM, DD, MJD
      integer t1901,cal(0:13,0:1),leap,j
      real*8  year
      parameter (t1901=15384,year=365.25d0)
      save cal
      data cal /-31,0,31,59,90,120,151,181,212,243,273,304,334,365,
     |          -31,0,31,60,91,121,152,182,213,244,274,305,335,366/
      j=yy
      if (j.lt.55) then
         j=j+100
      else if (j.gt.1900) then
         j=j-1900
      endif
      leap=0
      if (mod(j,4).eq.0) leap=1
      mjd=t1901+int((j-1)*year)+cal(mm,leap)+dd
      end

**MJD2YMD -- Convert MJD to YY, MM and DD
*+
      SUBROUTINE MJD2YMD (MJD, YY, MM, DD)
      implicit none
      INTEGER MJD, YY, MM, DD
      integer t1901,cal(0:13,0:1),leap,t
      real*8  year
      parameter (t1901=15384,year=365.25d0)
      save cal
      data cal /-31,0,31,59,90,120,151,181,212,243,273,304,334,365,
     |          -31,0,31,60,91,121,152,182,213,244,274,305,335,366/
      t=mjd-t1901-1
      yy=idint((t+8d-1)/year)
      t=t-int(yy*year)
      yy=yy+1901
      mm=t/31+1
      leap=0
      if (mod(yy,4).eq.0) leap=1
      if (t.ge.cal(mm+1,leap)) mm=mm+1
      dd=t-cal(mm,leap)+1
      end
