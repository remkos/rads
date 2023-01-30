/* Driver program for webtide subroutines

	 usage: webtide_c webtide_dir < input_file
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "webtide.h"

/* This defines the (Julian) date of the switch to the Gregorian calendar */
/*   (Used in the julday subroutine) */
#define IGREG (15 + 31L * (10 + 12L * 1582))

static int daytable[2][13] = {  /* Table of # days for each month */
	{0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
	{0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

long int julday (long int, long int, long int);
void get_date (int, int, long int *, long int *);

/*****************************************************************************/
int main (int argc, char *argv[])
{
	int nread;
	long int kd, year, month, day, hour, minute, dayofyear, cplx;
	double time, longitude, latitude, tide, tide2, seconds;
	struct WebTideInfo W;

	cplx = WebTideInit (argv[1], "s2c", &W);

/* Loop through the input file.
	 For each line calculate the tidal correction */
	while ((nread = fscanf (stdin, "%lf %lf %ld %ld %ld %ld %lf", &longitude, &latitude, &year, &dayofyear, &hour, &minute, &seconds)) != EOF) {
		if (nread != 7) {
			fprintf (stderr, "Error reading the input file!!\n");
			exit (1);
		}

		get_date (year, dayofyear, &month, &day);
		kd = julday (day, month, year) - 2446067;
		time = (double)kd * 86400.0 + (double)hour * 3600.0 + (double)minute * 60.0 + (double)seconds;
		WebTide (&W, &time, &latitude, &longitude, &tide, &tide2);

/* Ouput the tidal correction */
		if (cplx == 0) {
		fprintf (stdout, "%7.4lf %13.8lf %13.8lf %4ld %2ld %2ld %2ld %2ld %2ld %5.2lf\n", tide, longitude, latitude, year, month, day, dayofyear, hour, minute, seconds);
		}
		else {
		fprintf (stdout, "%7.4lf %7.4f %13.8lf %13.8lf %4ld %2ld %2ld %2ld %2ld %2ld %5.2lf\n", tide, tide2, longitude, latitude, year, month, day, dayofyear, hour, minute, seconds);
		}
	}

	WebTideFree (&W);

	return (0);
}

/*****************************************************************************/
void get_date (int year, int dayofyear, long int *month, long int *day)
/* Get the day and month from the day # of the year */
{
	int i, leap;

	leap =  (((year % 4 == 0) &&  (year % 100 != 0)) ||  (year % 400 == 0));
	for (i = 1; dayofyear > daytable[leap][i]; i++) {
		dayofyear -= daytable[leap][i];
	}
	*month = (long)i;
	*day = (long)dayofyear;
}

/*****************************************************************************/
long int julday (long int id, long int im, long int iy)
/* Calculate the Julian day number.
	 Accounts for the change to the Gregorian calandar. */
{
	long int jul, ja, jy, jm;

	jy = iy;
	if (jy == 0) {
		fprintf (stderr, "JULDAY: There is no year 0!\n");
		exit (-1);
	}
	if (jy < 0) ++jy;
	if (im > 2) {
		jm = im + 1;
	}
	else {
		--jy;
		jm = im + 13;
	}

	jul = (long int) (floor (365.25 * jy) + floor (30.6001 * jm) + id + 1720995);
	if ((id + 31L *  (im + 12L * iy)) >= IGREG) {
		ja = (long int) (0.01 * jy);
		jul += 2 - ja + (long int) (0.25 * ja);
	}
	return jul;
}
