      program mdate_main

* Convert almost any type of date into others.
* Examples of the input value:
* YYMMDD.DDD       = Year-month-day and fractions of day
* YYYYMMDD.DDD     = idem
* MJD.DDD          = Modified Julian Day and fractions of day
* YYDDD.DDD        = Year, day of the year and fractions of day
* YYYYDDD.DDD      = idem
* YYMMDDHHMMSS.SSS = Year-month-day hour-minute-seconds and fractions of
*                    seconds
* YYYYMMDDHHMMSS.S = idem
* SEC85.SSS        = UTC seconds since 1985 and fraction of seconds
*
* Requires:
* mdate, sec85, chrdat
*-
*   Jan-1996 - Created by Remko Scharroo
* 5-Aug-1998 - YYYYMMDD added, format flags added
*-----------------------------------------------------------------------
      implicit none
      character arg*80,cdate*15,ddate*15
      real*8    date/1d30/,sec,day,sec85,mjd85,mjd,ddd
      integer   yy,mdate,iargc,iarg,mode/0/
      parameter (day=86400d0,mjd85=46066d0)

      if (iargc().lt.1) goto 1300
      do iarg=1,iargc()
         call getarg(iarg,arg)
	 if (arg(1:2).eq.'-m') then
	    mode=1
	 else if (arg(1:2).eq.'-y') then
	    mode=2
	 else if (arg(1:2).eq.'-d') then
	    mode=3
	 else if (arg(1:2).eq.'-l') then
	    mode=4
	 else if (arg(1:2).eq.'-s') then
	    mode=5
	 else if (arg(1:1).eq.'-') then
	    goto 1300
	 else
	    read (arg,*) date
	 endif
      enddo
      if (date.gt.1d20) goto 1300
	    
      if (mode.eq.5) then
         sec=date
      else
      	 sec=sec85(mode,date)
      endif
      mjd=sec/day+mjd85
      call chrdat(nint(sec),cdate)
      call chrdat(int(sec),ddate)
      yy=mdate(1,int(mjd))/10000
      ddd=mjd-mdate(2,yy*10000+0101)+1
      write (*,600) sec,mjd,yy*1000+ddd,ddate(1:6),ddate(8:9),
     |	ddate(11:13),ddate(14:15),sec-int(sec),cdate
      goto 9999

1300  write (0,610)
600   format (
     |'          SEC85 =',f15.3/
     |'            MJD =',f13.6/
     |'          YYDDD =',f13.6/
     |'   YYMMDDHHMMSS = ',a6,3a2,f4.3/
     |'YYMMDD HH:MM:SS = ',a)
610   format ('mdate: convert almost any date format to others'//
     |'syntax: mdate [-format] <date>'//
     |'where <date> can be of the form:'/
     |' YYMMDD.DDD       = Year-month-day and fractions of day'/
     |' YYYYMMDD.DDD     = idem'/
     |' YYDDD.DDD        = Year, day of the year and fractions of day'/
     |' YYYYDDD.DDD      = idem'/
     |' MJD.DDD          = Modified Julian Day and fractions of day'/
     |' YYMMDDHHMMSS.SSS = year-month-day hour-minute-seconds and',
     |' fractions of seconds'/
     |' YYYYMMDDHHMMSS.S = idem'/
     |' SEC85.SSS        = UTC seconds since 1985 and fraction of',
     |' seconds'//
     |'and [-format] specifies the date format (if not specified, it',
     |' is guessed'/
     |' -y  = YYMMDD.DDD or YYYYMMDD.DD'/
     |' -l  = YYMMDDHHMMSS.S or YYYYMMDDHHMMSS.S'/
     |' -d  = YYDDD.DDD or YYYYDDD.DD'/
     |' -m  = MJD.DDD'/
     |' -s  = SEC85.SSS')
9999  end
