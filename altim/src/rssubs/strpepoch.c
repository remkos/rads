/*STRPEPOCH -- Fortran mimic of the C strptime routine
 +
      FUNCTION STRPEPOCH (BUF, FMT, EPOCH)
      FUNCTION STRP1970 (BUF, FMT)
      FUNCTION STRP1985 (BUF, FMT)
      FUNCTION STRP2000 (BUF, FMT)
      CHARACTER*(*) BUF, FMT
      INTEGER*4 EPOCH

  This function will make it easy on you to convert almost any time string
  into seconds from 1970, 1985, 2000 or any other epoch. The time string
  contained in BUF should have a format defined by FMT
  (see 'Format specifications' below).

  The function returns an integer value that counts the number of UTC
  seconds since 1 January 00:00 UTC of the year 1970, 1985 or 2000 or
  since the epoch defined by the value EPOCH (in seconds since 1 Jan 1970).
  Note that when only the date is given, time will be 00:00:00, when
  only time is given, the date will be 1 January or the given year.

  STRP* is basically the reverse of STRF*.

  Limitation: FMT and BUF will be truncated at 320 characters.

  Arguments:

      BUF     (input) : Time string that should conform to the format
                        specified in FMT
      FMT     (input) : Format specification of the time string
                        (see below)
      EPOCH   (input) : Epoch in seconds since 1 Jan 1970.
      STRPEPOCH (out) : Number of seconds since EPOCH.
      STRP1970  (out) : Number of seconds since 1 Jan 1970.
      STRP1985  (out) : Number of seconds since 1 Jan 1985.
      STRP2000  (out) : Number of seconds since 1 Jan 2000.

  Format specifications:

      A full list of the format specifications can be found in the manual
      of strftime ('man strftime'). This manual will give the exact list
      appropriate to your system. Nevertheless, we provide a list of the
      special sequences that can be used. But you can also use normal
      characters.

      %%  same as %
      %a  locale's abbreviated weekday name
      %A  locale's full weekday name
      %b  locale's abbreviated month name
      %B  locale's full month name
      %c  locale's appropriate date and time representation
      %C  locale's century number (the year divided by 100 and
          truncated to an integer) as a decimal number [00-99]
      %d  day of month ( 01 - 31 )
      %D  date as %m/%d/%y
      %e  day of month (1-31; single digits are preceded by a blank)
      %h  locale's abbreviated month name.
      %H  hour ( 00 - 23 )
      %I  hour ( 01 - 12 )
      %j  day number of year ( 001 - 366 )
      %KC locale's appropriate date and time representation
      %m  month number ( 01 - 12 )
      %M  minute ( 00 - 59 )
      %n  same as new-line
      %p  locale's equivalent of either AM or PM
      %r  time as %I:%M:%S [AM|PM]
      %R  time as %H:%M
      %S  seconds ( 00 - 61 ), allows for leap seconds
      %t  same as a tab
      %T  time as %H:%M:%S
      %U  week number of year ( 00 - 53 ), Sunday is the first day of week 1
      %w  weekday number ( 0 - 6 ), Sunday = 0
      %W  week number of year ( 00 - 53 ), Monday is the first day of week 1
      %x  locale's appropriate date representation
      %X  locale's appropriate time representation
      %y  year within century ( 00 - 99 )
      %Y  year as ccyy ( e.g. 1986)
      %Z  time zone name or no characters if no time zone exists

      Example: FMT="Today (%A %d/$b) I feel 100%%"
--
  12-Dec-2007 - Created by Remko Scharroo
------------------------------------------------------------------------
*/
/* Next two lines needed to make strptime and timegm work on GNU/Linux */
#define _XOPEN_SOURCE
#define _BSD_SOURCE

#include <string.h>
#include <time.h>
#include <sysdep.h>

#ifdef CAPITALS
#define strpepoch STRPEOCH
#define strp1970 STRP1970
#define strp1985 STRP1985
#define strp2000 STRP2000
#endif

#ifdef UNDERSCOREBEFORE
#define strpepoch _strpepoch
#define strp1970 _strp1970
#define strp1985 _strp1985
#define strp2000 _strp2000
#endif

#ifdef UNDERSCOREAFTER
#define strpepoch strpepoch_
#define strp1970 strp1970_
#define strp1985 strp1985_
#define strp2000 strp2000_
#endif

#define MIN(a,b) ((a) < (b) ? (a) : (b))

long strpepoch (char *buf, char *fmt, long *epoch, int bufsz, int fmtsz);

long strp1970(char *buf, char *fmt, int bufsz, int fmtsz)
{
	long epoch = 0;
	return strpepoch (buf, fmt, &epoch, bufsz, fmtsz);
}

long strp1985(char *buf, char *fmt, int bufsz, int fmtsz)
{
	long epoch = 473385600;
	return strpepoch (buf, fmt, &epoch, bufsz, fmtsz);
}

long strp2000(char *buf, char *fmt, int bufsz, int fmtsz)
{
	long epoch = 946684800;
	return strpepoch (buf, fmt, &epoch, bufsz, fmtsz);
}

long strpepoch(char *buf, char *fmt, long *epoch, int bufsz, int fmtsz)
{
	 char temp1[320], temp2[320];
	 struct tm time;
	 int l;
	 l = MIN(bufsz,320) - 1;
	 strncpy(temp1, buf, l+1);
	 while (temp1[l] == ' ' && l >= 0) l--;
	 temp1[l+1]='\0';
	 l = MIN(fmtsz,320) - 1;
	 strncpy(temp2, fmt, l+1);
	 while (temp2[l] == ' ' && l >= 0) l--;
	 temp2[l+1]='\0';
	 time = *gmtime(epoch);
	 strptime(temp1, temp2, &time);
	 return timegm(&time) - *epoch;
}
