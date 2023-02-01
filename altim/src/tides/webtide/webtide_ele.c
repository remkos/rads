/* Driver program for webtide subroutines */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "webtide.h"

int basis2d (double, double, double *);

/*****************************************************************************/
int main( int argc, char *argv[] )
{
  int i, j, k, cplx;
  double x, y, basis[3], t = 0.0, tide1, tide2;

  cplx = WebTideInit("webtide.cfg");

  for ( i = 0; i < 201; i++) {
     x = -66.00 + 0.01 * i;
     for ( j = 0; j < 201; j++) {
        y = 44.00 + 0.01 * j;
	k = basis2d (x, y, basis);
	k = WebTide (&t, &y, &x, &tide1, &tide2);
	fprintf (stdout, "%9.3f %9.3f %9.3f\n", x, y, tide1);
     }
  }

  WebTideFree ();

  return( 0 );
}
