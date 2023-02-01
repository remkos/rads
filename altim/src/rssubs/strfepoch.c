/*STRFEPOCH -- Fortran mimic of the C strftime routine
 +
      INTEGER*4 FUNCTION STRFEPOCH (BUF, FMT, SEC, EPOCH)
      INTEGER*4 FUNCTION STRF1970 (BUF, FMT, SEC)
      INTEGER*4 FUNCTION STRF1985 (BUF, FMT, SEC)
      INTEGER*4 FUNCTION STRF2000 (BUF, FMT, SEC)
      CHARACTER*(*) BUF, FMT
      INTEGER*4     SEC, EPOCH

  The function STRF1970 makes it easy on you to convert your system time
  (reckonned in seconds from 1 Jan 1970) into almost any sensible
  string representation (see 'Format specifications' below).
  STRF1985 and STRF2000 are also provided to do the same, but providing
  the amount of seconds since 1 Jan 1985 and 1 Jan 2000, respectively.
  Use STRFEPOCH to add any number of seconds first. EPOCH is the amount
  of seconds since 1 Jan 1970 to be added first.

  At the function call the integer variable SEC will contain the number
  of UTC seconds since the defined EPOCH (or 1970, 1985 or 2000). This
  time stamp is then transfered into a character string BUF based on
  the format specification FMT.

  The function returns an integer value that determines the lengths
  of the string after substitution according to the format, desregarding
  trailing spaces.

  If during substitution the capacity of the string buffer BUF is exceeded,
  the function returns a value 0 and the buffer will be empty (all spaces).

  STRP* is basically the reverse of STRF*.

  Limitation: FMT and BUF will be truncated at 320 characters.

  Arguments:

      BUF    (output) : Output string. The lengths of this string is
                        returned as the value of the function STRF1970
      FMT     (input) : Format specification of the output string
                        (see below)
      SEC     (input) : Clock value in seconds from 1970, 1985, or 2000
      EPOCH   (input) : Additional seconds since 1970.
      STRF*  (output) : Functions all return the number of characters of BUF.

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
  26-Feb-1999 - Created by Remko Scharroo
  24-Aug-1999 - Fix string overflow
  12-Dec-2007 - Changed from strf1970.c and strftime.F to strfepoch.c
   4-Mar-2008 - Fixed problem with Fortran wrapper on some systems:
                use int i.s.o. long
------------------------------------------------------------------------
*/
#include <string.h>
#include <time.h>
#include <sysdep.h>

#ifdef CAPITALS
#define strfepoch STRFEOCH
#define strf1970 STRF1970
#define strf1985 STRF1985
#define strf2000 STRF2000
#endif

#ifdef UNDERSCOREBEFORE
#define strfepoch _strfepoch
#define strf1970 _strf1970
#define strf1985 _strf1985
#define strf2000 _strf2000
#endif

#ifdef UNDERSCOREAFTER
#define strfepoch strfepoch_
#define strf1970 strf1970_
#define strf1985 strf1985_
#define strf2000 strf2000_
#endif

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int strfepoch(char *buf, char *fmt, int *sec, int *epoch, int bufsz, int fmtsz);

int strfepoch(char *buf, char *fmt, int *sec, int *epoch, int bufsz, int fmtsz)
{
	char temp1[320], temp2[320];
	int l, len;
	time_t clock;
	strncpy(temp1, fmt, fmtsz);
	for (l = MIN(fmtsz,320) - 1; temp1[l] == ' ' && l >= 0; l--) {};
	temp1[l+1]='\0';
	clock = (time_t)(*sec + *epoch);
	len = strftime(temp2, 320, temp1, gmtime(&clock));
	strncpy(buf, temp2, MIN(len,bufsz));
	for (l = len; l < bufsz; l++) buf[l]=' ';
	return len;
}

int strf1970(char *buf, char *fmt, int *sec, int bufsz, int fmtsz)
{
	int epoch = 0;
	return strfepoch(buf, fmt, sec, &epoch, bufsz, fmtsz);
}

int strf1985(char *buf, char *fmt, int *sec, int bufsz, int fmtsz)
{
	int epoch = 473385600;
	return strfepoch(buf, fmt, sec, &epoch, bufsz, fmtsz);
}

int strf2000(char *buf, char *fmt, int *sec, int bufsz, int fmtsz)
{
	int epoch = 946684800;
	return strfepoch(buf, fmt, sec, &epoch, bufsz, fmtsz);
}
