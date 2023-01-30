/*FDATE -- Get date and time as character string
 +
      SUBROUTINE FDATE (STRING)
      CHARACTER*(*) STRING
 
  Return the current date and time.
  To receive the whole string, the STRING should be declared
  at least CHARACTER*24 in the calling (sub)program.
 
  Argument:
   STRING : receives date and time, truncated or extended with
            blanks as necessary.
 -
  5-Jul-1993 - New manual
  5-Jan-1996 - Multi-platform version
 -
*/

#include <sys/types.h>
#include <time.h>
#include <string.h>

#include <sysdep.h>

#ifdef CAPITALS
#define fdate FDATE
#endif
#ifdef UNDERSCOREBEFORE
#define fdate _fdate
#endif
#ifdef UNDERSCOREAFTER
#define fdate fdate_
#endif

void fdate(string, string_len)
char *string;
long string_len;
{
   long now,
      i;
   char *time_string;

   now = time((time_t *) 0);
   time_string = ctime(&now);

   if (string_len <= 24L) {
      memcpy(string, time_string, string_len);
   } else {
      memcpy(string, time_string, 24);
      for (i = 24L; i < string_len; i++) {
         string[i] = ' ';
      }
   }
}
