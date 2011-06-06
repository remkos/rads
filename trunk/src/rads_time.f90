module rads_time

contains

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

!*sec85ymdhms -- Convert SEC85 to YY, MM, DD, HH, MN, SS
!+
pure subroutine sec85ymdhms (sec, yy, mm, dd, hh, mn, ss)
use typesizes
real(eightbytereal), intent(in) :: sec
integer(fourbyteint), intent(out) :: yy, mm, dd, hh, mn
real(eightbytereal), intent(out) :: ss
!
! Converts SEC85 (seconds since 1985-01-01 00:00:00) to year, month,
! day, hours, minutes, and (floating point) seconds.
!
! Arguments:
!  sec : Seconds since 1 Jan 1985
!  yy  : 4-digit year (between 1901 and 2099)
!  mm  : Month
!  dd  : Day
!  hh  : Hours
!  mn  : Minutes
!  ss  : Seconds (floating point)
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
end subroutine sec85ymdhms

!*strf1985f -- Construct date string from seconds since 1985
!+
elemental function strf1985f (sec, sep)
use typesizes
real(eightbytereal), intent(in) :: sec
character(len=1), optional, intent(in) :: sep
character(len=26) :: strf1985f
!
! This routine formats relative time in seconds since 1985 into a
! character string of 26 bytes in the form YYYY-MM-DDxHH:MM:SS.SSSSSS,
! where 'x' is a seperator character given by the optional argument
! <sep> (default is a space)
!
! Arguments:
!  sec       : Seconds since 1.0 Jan 1985
!  sep       : Seperator character (default is a space)
!
! Return value:
!  strf1985f : Character string of time
!-
! Copyright (c) Remko Scharroo, Altimetrics LLC
!-----------------------------------------------------------------------
integer(fourbyteint) :: yy, mm, dd, hh, mn
real(eightbytereal) :: ss
character(len=1) :: x
call sec85ymdhms (sec, yy, mm, dd, hh, mn, ss)
if (present(sep)) then
	x = sep
else
	x = ' '
endif
write (strf1985f, '(i4.4,2("-",i2.2),a1,2(i2.2,":"),i2.2,f7.6)') &
	yy, mm, dd, x, hh, mn, floor(ss), modulo(ss,1d0)
end function strf1985f

!*strp1985f -- Parse date string and convert it to seconds since 1985
!+
elemental function strp1985f (string)
use typesizes
character(len=26), intent(in) :: string
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
integer(fourbyteint) :: yy,mm,dd,hh,mn,mjd,ios
real(eightbytereal) :: ss
hh = 0 ; mn = 0 ; ss = 0d0
read (string,'(i4,4(1x,i2),1x,f9.6)',iostat=ios) yy,mm,dd,hh,mn,ss
call ymd2mjd(yy,mm,dd,mjd)
strp1985f = (mjd - 46066) * 86400d0 + hh * 3600d0 + mn * 60d0 + ss
end function strp1985f

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
! The first parameter (i) defines the input format: i=1 indicates MJD,
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
!-
! 1995 -- Created by Remko Scharroo
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
	ff = date - dd*1000000
	yy = dd/10000
	dd = dd - yy*10000
	mm = dd/100
	dd = dd - mm*100
	hh = ff/1d4
	ff = ff - hh*1d4
	mn = ff/1d2
	ff = ff - mn*1d2
	call ymd2mjd (yy,mm,dd,mjd)
	sec85 = (mjd-mjd85)*day + hh*3600 + mn*60 + ff
case default
	sec85 = date
end select

end function sec85

!*datearg -- Processes standard datation arguments
!+
function datearg (arg, t0, t1, dt)
use typesizes
logical :: datearg
character*(*), intent(in) :: arg
real(eightbytereal), intent(out) :: t0
real(eightbytereal), intent(out), optional :: t1, dt
!
! This function processes standard datation arguments of either of the
! following forms:
!
!   mjd=t0[,t1][,dt] : Modified Julian Dates
!   sec=t0[,t1][,dt] : Seconds since 1.0 Jan 1985
!   ymd=t0[,t1][,dt] : [YY]YYMMDD.DDD or [YY]YYMMDDHHMMSS.SSS
!   doy=t0[,t1][,dt] : YYDDD.DDD
!     t=t0[,t1][,dt] : MJD.DDD or [YY]YYMMDD.DDD or [YY]YYMMDDHHMMSS.SSS
!
! Normally <arg> is read from the argument line of a command using the getarg
! routine. When <arg> is parsed through to datearg it checks whether <arg> is
! of one of the above forms. If not, datearg gets the value .false., and
! none of the arguments is altered.
!
! If <arg> is of one of the standard datation forms, the values of <t0> and
! <t1> (when given) are read and interpreted as stated above and converted
! to seconds since 1.0 Jan 1985. When provided, <dt> is read too, but is not
! converted. On return, datearg will get the value .true.
!
! When <t1> or <dt> are not provided in <arg>, the arguments of the function
! call will keep their values unchanged.
! Since <t1> and <dt> are optional arguments, they can also be omitted from
! the call to datearg.
!
! Arguments:
!   datearg (output): .true. if <arg> is in a standard datation format
!   arg      (input): datation argument (see above)
!   t0, t1  (output): datation arguments converted to SEC85
!                     (seconds since 1.0 Jan 1985)
!   dt      (output): third datation argument (not converted)
!-
! 03-Nov-1999 : Created by Remko Scharroo
!-----------------------------------------------------------------------
integer :: ios,l,mode
real(eightbytereal) :: tt0,tt1

! Scan the argument for datation format

l = index(arg,'=')
if (l /= 2 .and. l /= 4) then
	datearg = .false.
	return
else if (arg(:2) == 't=' .or. arg(:2) == 'T=') then
	mode=0
else if (arg(:4) == 'sec=' .or. arg(:4) == 'SEC=') then
	mode=-1
else if (arg(:4) == 'mjd=' .or. arg(:4) == 'MJD=') then
	mode=1
else if (arg(:4) == 'doy=' .or. arg(:4) == 'DOY=') then
	mode=3
else if (arg(:4) == 'ymd=' .or. arg(:4) == 'YMD=') then
	mode=5
else
	datearg = .false.
	return
endif

! Initialize the values to an unlikely value

tt0 = 1d50
tt1 = 1d50

! Read the values, and convert when specified

if (present(t1) .and. present(dt)) then
	read (arg(l+1:),*,iostat=ios) tt0,tt1,dt
	if (abs(tt1) < 1d49) t1 = sec85(mode,tt1)
else if (present(t1)) then
	read (arg(l+1:),*,iostat=ios) tt0,tt1
	if (abs(tt1) < 1d49) t1 = sec85(mode,tt1)
else if (present(dt)) then
	read (arg(l+1:),*,iostat=ios) tt0,dt
else
	read (arg(l+1:),*,iostat=ios) tt0
endif
if (abs(tt0) < 1d49) t0 = sec85(mode,tt0)
datearg = .true.
end function datearg

end module rads_time
