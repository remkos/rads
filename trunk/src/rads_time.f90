!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2014  Remko Scharroo (Altimetrics LLC)
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!-----------------------------------------------------------------------

module rads_time

!***********************************************************************
!*sec85ymdhms -- Convert SEC85 to YY, MM, DD, HH, MN, SS
!+
! pure subroutine sec85ymdhms (sec, yy, mm, dd, hh, mn, ss)
! use typesizes
! real(eightbytereal) <or> integer(fourbyteint), intent(in) :: sec
! integer(fourbyteint), intent(out) :: yy, mm, dd, hh, mn
! real(eightbytereal) <or> integer(fourbyteint), intent(out) :: ss
!
! Converts SEC85 (seconds since 1985-01-01 00:00:00) to year, month,
! day, hours, minutes, and seconds.
!
! The input <sec> can be either a floating point value (double) or an
! integer. The return value for <ss> will be likewise a double or an
! integer.
!
! Arguments:
!  sec : Seconds since 1 Jan 1985 (double float or integer)
!  yy  : 4-digit year (between 1901 and 2099)
!  mm  : Month
!  dd  : Day
!  hh  : Hours
!  mn  : Minutes
!  ss  : Seconds (double float or integer)
!-----------------------------------------------------------------------
private :: sec85ymdhms_dble, sec85ymdhms_int4
interface sec85ymdhms
	module procedure sec85ymdhms_dble
	module procedure sec85ymdhms_int4
end interface sec85ymdhms

!***********************************************************************
!*strf1985f -- Construct date string from seconds since 1985
!+
! elemental function strf1985f (sec, sep, frac)
! use typesizes
! real(eightbytereal) <or> integer(fourbyteint), intent(in) :: sec
! character(len=1), optional, intent(in) :: sep
! logical, optional, intent(in) :: frac
! character(len=26) <or> character(len=19) :: strf1985f
!
! This routine formats relative time in seconds since 1985 into a
! character string of either:
! 1) 26 bytes in the form YYYY-MM-DDxHH:MM:SS.SSSSSS (when <sec> is a float)
! 2) 19 bytes in the form YYYY-MM-DDxHH:MM:SS (when <sec> is an integer)
! where 'x' is a separator character given by the optional argument
! <sep> (default is a space)
!
! Arguments:
!  sec       : Seconds since 1.0 Jan 1985 (double float or integer)
!  sep       : (Optional) Separator character (default is a space)
!
! Return value:
!  strf1985f : Character string of time
!-----------------------------------------------------------------------
private :: strf1985f_dble, strf1985f_int4
interface strf1985f
	module procedure strf1985f_dble
	module procedure strf1985f_int4
end interface strf1985f

contains

!***********************************************************************
!*ymd2mjd -- Convert YY, MM and DD to MJD
!+
pure subroutine ymd2mjd (yy, mm, dd, mjd)
use typesizes
integer(fourbyteint), intent(in) :: yy, mm, dd
integer(fourbyteint), intent(out) :: mjd
!
! Converts year, month and date to modified julian dates.
! This routine works both with year numbers in one of two forms:
! 2-digit numbers: assumes the year is between 1955 and 2054.
! 4-digit years: either 19xx or 20xx (works for 1901-2099)
!
! For efficiency, month numbers 0 and 13 are allowed (for example
! when you want to increase or decrease a date by one month).
! No checks are made on the validity of month or day numbers.
!
! This routine can also be used to convert year and day-of-year to
! MJD by setting the month number to 1.
!
! Arguments:
!  yy  : Year (either as 2 digits or 4 digits)
!  mm  : Month
!  dd  : Day (or day-of-year)
!  mjd : Modified Julian Day
!-----------------------------------------------------------------------
integer(fourbyteint), parameter :: mjd1901 = 15385 ! 1 Jan 1901
integer(fourbyteint), parameter :: cal(0:13) = &
	(/-31,0,31,59,90,120,151,181,212,243,273,304,334,365/)
integer (fourbyteint) :: j
! j is number of years since 1901
if (yy < 55) then
	j = yy + 99
else if (yy > 1900) then
	j = yy - 1901
else
	j = yy - 1
endif
mjd = mjd1901 + j*365 + j/4 + cal(mm) + (dd - 1)
if (modulo(yy,4) == 0 .and. mm > 2) mjd = mjd + 1
end subroutine ymd2mjd

!*mjd2ymd -- Convert MJD to YY, MM and DD
!+
pure subroutine mjd2ymd (mjd, yy, mm, dd)
use typesizes
integer(fourbyteint), intent(in) :: mjd
integer(fourbyteint), intent(out) :: yy, mm, dd
!
! Converts Modified Julian Date to year, month and day.
!
! This routine works only for the years 1901 through 2099.
!
! Arguments:
!  mjd : Modified Julian Day
!  yy  : 4-digit year (between 1901 and 2099)
!  mm  : Month
!  dd  : Day (or day-of-year when no <mm> is given)
!-----------------------------------------------------------------------
integer (fourbyteint) :: leap
integer(fourbyteint), parameter :: cal(0:13,0:1) = &
	reshape((/-31,0,31,59,90,120,151,181,212,243,273,304,334,365, &
	          -31,0,31,60,91,121,152,182,213,244,274,305,335,366/),(/14,2/))
call mjd2yd (mjd, yy, dd)
if (modulo(yy,4) == 0) then
	leap = 1
else
	leap = 0
endif
mm = dd/32 + 1
if (dd > cal(mm+1,leap)) mm = mm + 1
dd = dd - cal(mm,leap)
end subroutine mjd2ymd

!***********************************************************************
!*mjd2yd -- Convert MJD to YY and DDD
!+
pure subroutine mjd2yd (mjd, yy, ddd)
use typesizes
integer(fourbyteint), intent(in) :: mjd
integer(fourbyteint), intent(out) :: yy, ddd
!
! Converts Modified Julian Date to year and day-of-year.
!
! This routine works only for the years 1901 through 2099.
!
! Arguments:
!  mjd : Modified Julian Day
!  yy  : 4-digit year (between 1901 and 2099)
!  ddd : Day-of-year
!-----------------------------------------------------------------------
integer (fourbyteint) :: j
integer(fourbyteint), parameter :: mjd1901 = 15385 ! 1 Jan 1901
ddd = mjd - mjd1901 ! Number of days since 1 Jan 1901
j = floor((ddd+0.8d0)/365.25d0) ! Number of years since 1901
yy = j + 1901
ddd = ddd - j*365 - j/4 + 1 ! Day of year
end subroutine mjd2yd

!*sec85ymdhms_dble -- Convert SEC85 to YY, MM, DD, HH, MN, SS
!+
pure subroutine sec85ymdhms_dble (sec, yy, mm, dd, hh, mn, ss)
use typesizes
real(eightbytereal), intent(in) :: sec
integer(fourbyteint), intent(out) :: yy, mm, dd, hh, mn
real(eightbytereal), intent(out) :: ss
!
! Converts SEC85 (seconds since 1985-01-01 00:00:00) to year, month,
! day, hours, minutes, and (floating point) seconds.
!-----------------------------------------------------------------------
real(eightbytereal) :: tt
! First split off the day part
dd = floor(sec/86400d0)
tt = sec - dd * 86400
call mjd2ymd (dd + 46066, yy, mm, dd)
! Now handle the fractions of day
hh = floor(tt/3600d0)
tt = tt - hh * 3600
mn = floor(tt/60d0)
ss = tt - mn * 60
end subroutine sec85ymdhms_dble

!***********************************************************************
!*sec85ymdhms_int -- Convert SEC85 to YY, MM, DD, HH, MN, SS
!+
pure subroutine sec85ymdhms_int4 (sec, yy, mm, dd, hh, mn, ss)
use typesizes
integer(fourbyteint), intent(in) :: sec
integer(fourbyteint), intent(out) :: yy, mm, dd, hh, mn, ss
!
! Converts SEC85 (seconds since 1985-01-01 00:00:00) to year, month,
! day, hours, minutes, and (floating point) seconds.
!-----------------------------------------------------------------------
integer(fourbyteint) :: tt
! First split off the day part
dd = sec / 86400
tt = sec - dd * 86400
call mjd2ymd (dd + 46066, yy, mm, dd)
! Now handle the fractions of day
hh = tt / 3600
tt = tt - hh * 3600
mn = tt / 60
ss = tt - mn * 60
end subroutine sec85ymdhms_int4

!***********************************************************************
!*strf1985f_dble -- Construct date string from seconds since 1985
!+
elemental function strf1985f_dble (sec, sep)
use typesizes
real(eightbytereal), intent(in) :: sec
character(len=1), optional, intent(in) :: sep
character(len=26) :: strf1985f_dble
!-----------------------------------------------------------------------
integer(fourbyteint) :: yy, mm, dd, hh, mn
real(eightbytereal) :: ss
character(len=1) :: x
if (present(sep)) then
	x = sep
else
	x = ' '
endif
call sec85ymdhms (sec, yy, mm, dd, hh, mn, ss)
write (strf1985f_dble, '(i4.4,2("-",i2.2),a1,2(i2.2,":"),i2.2,f7.6)') &
		yy, mm, dd, x, hh, mn, floor(ss), modulo(ss,1d0)
end function strf1985f_dble

!***********************************************************************
!*strf1985f_int4 -- Construct date string from seconds since 1985
!+
elemental function strf1985f_int4 (sec, sep)
use typesizes
integer(fourbyteint), intent(in) :: sec
character(len=1), optional, intent(in) :: sep
character(len=19) :: strf1985f_int4
!-----------------------------------------------------------------------
integer(fourbyteint) :: yy, mm, dd, hh, mn, ss
character(len=1) :: x
if (present(sep)) then
	x = sep
else
	x = ' '
endif
call sec85ymdhms (sec, yy, mm, dd, hh, mn, ss)
write (strf1985f_int4, '(i4.4,2("-",i2.2),a1,2(i2.2,":"),i2.2)') &
		yy, mm, dd, x, hh, mn, ss
end function strf1985f_int4

!***********************************************************************
!*strp1985f -- Parse date string and convert it to seconds since 1985
!+
elemental function strp1985f (string)
use typesizes
character(len=*), intent(in) :: string
real(eightbytereal) :: strp1985f
!
! This routine reads a string of the form YYYY-MM-DDxHH:MM:SS.SSSSSS,
! where x is any character, and converts it to a seconds since
! 1.0 Jan 1985. Fractional seconds can be included or the HH:MM:SS
! part can be omitted entirely (to produce 00:00:00)
!
! Arguments:
!  string    : Character string of time
!
! Return value:
!  strp1985f : Seconds since 1.0 Jan 1985
!-----------------------------------------------------------------------
integer(fourbyteint) :: yy,mm,dd,hh,mn,ss,mjd,ios
real(eightbytereal) :: fs
hh = 0 ; mn = 0 ; ss = 0 ; fs = 0d0
read (string,'(i4,5(1x,i2))',iostat=ios) yy,mm,dd,hh,mn,ss
if (len(string) > 19) read (string(20:),*,iostat=ios) fs
call ymd2mjd(yy,mm,dd,mjd)
strp1985f = (mjd - 46066) * 86400d0 + hh * 3600d0 + mn * 60d0 + ss + fs
end function strp1985f

!***********************************************************************
!*sec85 -- Convert MJD or YYMMDD or YYDDD to SEC85
!+
pure function sec85 (i, date)
use typesizes
integer(fourbyteint), intent(in) :: i
real(eightbytereal), intent(in) :: date
real(eightbytereal) :: sec85
!
! This function converts Modified Julian Dates (MJD) or Year-Month-Day
! (YYMMDD) or Year-Day (YYDDD) or Year-Month-Day-Hours-Minutes-Seconds
! (YYMMDDHHMMSS) to seconds from 1.0 Jan 1985 (SEC85). Fractions of days or
! seconds can also be included.
!
! Except for a 2-digit year indication (YY) it is now also possible to use
! a 4-digit indication (YYYY). The routine automatically recognises the
! form: whether YY or YYYY is specified.
!
! The first parameter <i> defines the input format: i=1 indicates MJD,
! i=2 means YYMMDD or YYYYMMDD, i=3 specifies YYDDD or YYYYDDD, and
! i=4 implies YYMMDDHHMMSS or YYYYMMDDHHMMSS, i=5 is a combination of i=2
! and i=4 (i.e, it takes YYMMDD, YYYYMMDD, YYMMDDHHMMSS and YYYYMMDDHHMMSS,
! as input).
!
! One can also choose NOT to specify the input format and let the function
! 'guess' which format it is. Obviously, there are limitations. it recognises
!         YYMMDD :  000101 - 040101  ( 2000/01/01 - 2004/01/01 )
!            MJD :   40587 - 53005   ( 1970/01/01 - 2004/01/01 )
!          YYDDD :   53005 - 99365   ( 1953/01/05 - 2000/01/01 )
!         YYMMDD :  500101 - 991232  ( 1950/01/01 - 2000/01/01 )
!        YYYYDDD :  1950d3 - 2050d3  ( 1950/01/01 - 2050/01/01 )
!       YYYYMMDD :  1950d4 - 2050d4  ( 1950/01/01 - 2050/01/01 )
!   YYMMDDHHMMSS :     1d8 - 1d12    ( 1950/01/01 - 2050/01/01 )
! YYYYMMDDHHMMSS : 1950d10 - 2050d10 ( 1950/01/01 - 2050/01/01 )
!
! Arguments:
!  i     : i=1, convert MJD to SEC85
!          i=2, convert YYMMDD or YYYYMMDD to SEC85
!          i=3, convert YYDDD or YYYYDDD to SEC85
!          i=4, convert YYMMDDHHMMSS or YYYYMMDDHHMMSS to SEC85
!          i=5, convert [YY]YYMMDD[HHMMSS] to SEC85
!          i=0, convert any of the above to SEC85
!  date  : Input value.
!  sec85 : Output value, seconds from 1.0 Jan 1985.
!-----------------------------------------------------------------------
integer(fourbyteint) :: mjd,yy,mm,dd,hh,mn,ii
real(eightbytereal) :: ff
integer, parameter :: mjd85 = 46066, day = 86400

select case (i)
case (5) ! YYMMDD or YYYYMMDD or YYMMDDHHMMSS or YYYYMMDDHHMMSS -> SEC85
	if (date < 1d8) then
		ii = 2
	else
		ii = 4
	endif
case (0) ! Guess
	if (date < 0) then
		sec85 = date
		return
	else if (date < 40587d0) then
		ii = 2
	else if (date < 53005d0) then
		ii = 1
	else if (date < 500101d0) then
		ii = 3
	else if (date < 1950d3) then
		ii = 2
	else if (date < 1950d4) then
		ii = 3
	else if (date < 1d8) then
		ii = 2
	else if (date < 2050d10) then
		ii = 4
	else
		sec85 = date
		return
	endif
case default
	ii = i
end select

select case (ii)
case (1) ! MJD -> SEC85
	sec85 = (date-mjd85)*day
case (2) ! YYMMDD or YYYYMMDD -> SEC85
	dd = int(date)
	ff = date - dd
	yy = dd/10000
	dd = dd - yy*10000
	mm = dd/100
	dd = dd - mm*100
	call ymd2mjd (yy,mm,dd,mjd)
	sec85 = (mjd-mjd85+ff)*day
case (3) ! YYDDD or YYYYDDD -> SEC85
	yy = int(date/1d3)
	ff = date - yy*1000 - 1
	call ymd2mjd (yy,1,1,mjd)
	sec85 = (mjd-mjd85+ff)*day
case (4) ! YYMMDDHHMMSS or YYYYMMDDHHMMSS -> SEC85
	dd = int(date/1d6)
	ff = date - dd*1d6
	yy = dd/10000
	dd = dd - yy*10000
	mm = dd/100
	dd = dd - mm*100
	hh = floor(ff/1d4)
	ff = ff - hh*1d4
	mn = floor(ff/1d2)
	ff = ff - mn*1d2
	call ymd2mjd (yy,mm,dd,mjd)
	sec85 = (mjd-mjd85)*day + hh*3600 + mn*60 + ff
case default
	sec85 = date
end select

end function sec85

!***********************************************************************
!*dateopt -- Processes standard datation options
!+
function dateopt (optopt, optarg, t0, t1, dt)
use typesizes
logical :: dateopt
character(len=*), intent(in) :: optopt, optarg
real(eightbytereal), intent(out) :: t0
real(eightbytereal), intent(out), optional :: t1, dt
!
! This function processes standard datation arguments of either of the
! following forms:
!
!   --mjd=t0[,t1][,dt] : Modified Julian Dates
!   --sec=t0[,t1][,dt] : Seconds since 1.0 Jan 1985
!   --ymd=t0[,t1][,dt] : [YY]YYMMDD.DDD or [YY]YYMMDDHHMMSS.SSS
!   --doy=t0[,t1][,dt] : YYDDD.DDD
!  --time=t0[,t1][,dt] : MJD.DDD or [YY]YYMMDD.DDD or [YY]YYMMDDHHMMSS.SSS
!
! At input, these arguments have already been split up in an option part
! (without the --) and an argument part (after the =). This can be the
! result of the <getopt> routine, for example.
!
! When <optopt> is parsed through to <dateopt> it checks whether <optopt> is
! of one of the above forms. If not, <dateopt> gets the value .false., and
! none of the arguments is altered.
!
! If <optopt> is of one of the standard datation forms, the values of <t0> and
! <t1> (when given) are read from <optarg> and interpreted as stated above
! and converted to seconds since 1.0 Jan 1985. When provided, <dt> is read
! too, but is not converted. On return, <dateopt> will get the value .true.
! Values can be separated by spaced or slashes.
!
! When <t1> or <dt> are not provided in <optarg>, the arguments of the function
! call will keep their values unchanged.
! Since <t1> and <dt> are optional arguments, they can also be omitted from
! the call to <dateopt>.
!
! Arguments:
!   dateopt (output): .true. if <optopt> is in a standard datation format
!   optopt   (input): datation option ('mjd', 'sec', etc.) (see above)
!   optarg   (input): datation option argument (see above)
!   t0, t1  (output): datation arguments converted to SEC85
!                     (seconds since 1.0 Jan 1985)
!   dt      (output): third datation argument (not converted)
!-----------------------------------------------------------------------
integer :: i,mode
real(eightbytereal) :: tt0,tt1
character(len=len(optarg)) :: arg
real(eightbytereal), parameter :: nan = transfer ((/not(0_fourbyteint),not(0_fourbyteint)/),0d0)

! Check the options
select case (optopt)
case ('time', 't:')
	mode=0
case ('sec')
	mode=-1
case ('mjd')
	mode=1
case ('doy')
	mode=3
case ('ymd')
	mode=5
case default
	dateopt = .false.
	return
end select

! Initialize the values to NaN
tt0 = nan
tt1 = nan

! Replace any slashes with commas
arg = optarg
do i = 1,len_trim(optarg)
	if (arg(i:i) == '/') arg(i:i) = ','
enddo

! Read the values, and convert when specified
if (present(t1) .and. present(dt)) then
	read (arg, *, iostat=i) tt0, tt1, dt
	if (tt1 == tt1) t1 = sec85 (mode, tt1)
else if (present(t1)) then
	read (arg, *, iostat=i) tt0, tt1
	if (tt1 == tt1) t1 = sec85 (mode, tt1)
else if (present(dt)) then
	read (arg, *, iostat=i) tt0, dt
else
	read (arg, *, iostat=i) tt0
endif
if (tt0 == tt0) t0 = sec85 (mode, tt0)

! Successful return
dateopt = .true.
end function dateopt

!***********************************************************************
!*datestamp -- Create character string with current date
!+
function datestamp ()
character(len=10) :: datestamp
!
! This function produces the current date (in UTC) in the form
! 2012-06-15
!-----------------------------------------------------------------------
integer :: values(9)
call gmtime(time(),values)
write (datestamp, '(i4.4,2("-",i2.2))') values(6)+1900,values(5)+1,values(4)
end function datestamp

!***********************************************************************
!*timestamp -- Create character string with current date
!+
function timestamp (sep)
character(len=1), optional :: sep
character(len=19) :: timestamp
!
! This function produces the current date (in UTC) in the form
! 2012-06-15x02:58:15, where x is a separator character indicated
! by the optional argument <sep>. Default is a space.
!-----------------------------------------------------------------------
integer :: values(9)
character(len=1) :: x
if (present(sep)) then
	x = sep
else
	x = ' '
endif
call gmtime(time(),values)
write (timestamp, '(i4.4,2("-",i2.2),a1,i2.2,2(":",i2.2))') &
	values(6)+1900,values(5)+1,values(4),x,values(3:1:-1)
end function timestamp

end module rads_time
