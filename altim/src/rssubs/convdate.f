**CONVDATE -- Convert between time formats
*+
      program convdate

* This program converts time coordinates. Regularly, people will use one of
* the following conventions for denoting time:
*
* MJD              = Modified Julian Days
* YYDDD            = Year and day of year
* YYMMDD           = Year-Month-Day
* YYMMDDHHMMSS     = Year-Month-Day-Hours-Minutes-Seconds
* SEC85            = Seconds from 1.0 Jan 1985 (altimeter default)
*
* In order to facilitate this, the RSSUBS library contains several
* subroutines: MDATE, SEC85, CHRDAT and STRF1985. This program is an
* interface to these subroutines.
*
* Usage:
* convdate <date-specification> <format-specification>
* (run convdate without arguments for a full description)
*-
* 03-Nov-1999 - Created by Remko Scharroo
* 22-Feb-2000 - Option off?= added
* 03-Apr-2002 - date-int(date) > .9990 problem fixed for ymd output [ED/RS]
*-----------------------------------------------------------------------
      real*8 date/1d40/,sec85,dum,offset/0d0/,day,rel/0d0/
      integer*4 l,strf1985,iargc,iarg,mjd85,idate,msec
      character*80 arg,string,format/'%c'/
      logical*4 datearg
      parameter (day=86400d0,mjd85=46066)

* Scan arguments

      do iarg=1,iargc()
         call getarg(iarg,arg)
	 if (arg(:2).eq.'-h') then
	    goto 1300
	 else if (arg(:3).eq.'off' .and. arg(5:5).eq.'=') then
	    read (arg(6:),*) offset
	    if (arg(4:4).eq.'m') offset=offset*60d0
	    if (arg(4:4).eq.'h') offset=offset*3600d0
	    if (arg(4:4).eq.'d') offset=offset*day
	 else if (datearg(arg,date,dum,dum)) then
	 else if (arg(:3).eq.'rel' .and. datearg(arg(4:),rel,dum,dum))
     |		then
	 else if (arg.eq.'+YMD') then
	    format='%y%m%d%H%M%S'
	 else if (arg.eq.'+doy') then
	    format='%y%j'
	 else if (arg.eq.'+txt') then
	    format='%c'
	 else if (arg(:1).eq.'+') then
	    format=arg(2:)
	 else
	    read (arg,*) date
	    date=sec85(0,date)
	 endif
      enddo

* Add offset to date (now in SEC85)

      date=date+offset
      if (date.gt.1d30) goto 1300

* Print out date in requested format

      if (format.eq.'sec') then
         write (*,"(f14.3)") date-rel
      else if (format.eq.'SEC') then
         write (*,"(i10)") nint(date-rel)
      else if (format.eq.'day') then
         write (*,"(f14.3)") (date-rel)/day
      else if (format.eq.'DAY') then
         write (*,"(i10)") nint((date-rel)/day)
      else if (format.eq.'mjd') then
         write (*,"(f14.8)") date/day+mjd85
      else if (format.eq.'MJD') then
         write (*,"(i5)") nint(date/day)+mjd85
      else if (format.eq.'ymd') then
	 idate=nint(date-0.5d0)
	 msec=nint((date-idate)*1d3)
	 if (msec.gt.999) then
	    msec=msec-1000
	    idate=idate+1
	 endif
         l=strf1985(string,'%y%m%d%H%M%S',idate)
 	 write (*,"(a,'.',i3.3)") string(:l),msec
      else
         l=strf1985(string,format,nint(date))
	 write (*,"(a)") string(:l)
      endif

      goto 9999

1300  write (0,1310)
1310  format ('convdate - convert between time formats'//
     |'usage: convdate [options] <date> <format>'//
     |'Converts specified <date> into a string according to <format>'//
     |'where [options] are:'/
     |'  offs=delta : add ''delta'' seconds to given <date>'/
     |'  offm=delta : add ''delta'' minutes to given <date>'/
     |'  offh=delta : add ''delta'' hours to given <date>'/
     |'  offd=delta : add ''delta'' days to given <date>'/
     |'  rel<date>  : determine days or seconds relative to <date>',
     |' [Default: ymd=19850101]'//
     |'and <date> is either:'/
     |'  ymd=[YY]YYMMDD.D       : Year-Month-Day'/
     |'  ymd=[YY]YYMMDDHHMMSS.S : Year-Month-Day',
     |'-Hour-Minutes-Seconds'/
     |'  mjd=MJD.DDD            : Modified Julian Day'/
     |'  doy=[YY]YYDDD.DDD      : Year-DayOfYear'/
     |'  sec=SEC85.SSS          : Seconds since 1.0 Jan 1985'/
     |'  [t=]MJD or [YY]YYMMDD[HHMMSS] : free choise'//
     |'and <format> is either:'/
     |'  +ymd    : YYMMDDHHMMSS.SSS (real)'/
     |'  +YMD    : same as +%y%m%d%H%M%S (int)'/
     |'  +mjd    : Modified Julian Day (real)'/
     |'  +MJD    : Modified Julian Day (int)'/
     |'  +day    : Days since epoch (real)'/
     |'  +DAY    : Days since epoch (nint)'/
     |'  +sec    : Seconds since epoch (real)'/
     |'  +SEC    : Seconds since epoch (nint)'/
     |'  +txt    : same as +%c (int, default)'/
     |'  +doy    : same as +%y%j (int)'/
     |'  +<spec> : Format spec like ''date'' command (int)')

9999  end
