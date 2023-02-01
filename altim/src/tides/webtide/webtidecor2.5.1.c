/***********************************************************************
Copyright (c) Fisheries and Oceans Canada, 2010

This program may be freely redistributed for non-commercial purposes
under the condition that the copyright notices are not removed. You may
distribute modified versions of this program under the conditions that
both source and object code are made available without charge and that
clear notice is given of the modifications. Modified programs and source
code are not to be represented as being endorsed by, or made in
cooperation with, Fisheries and Oceans Canada.

This program is distributed by the copyright holder and contributors in
the hope that it will be useful. It is provided "AS IS" without
warranties of any kind, whether express or implied. This includes, but
is not limited to, merchantability or fitness for a particular purpose.

In no event shall the copyright holder be liable for damages, including
any general, special, incidental or consequential damages arising out of
the use or inability to use the program (including, but not limited to,
loss of use, data or profits).

This program is not to be used for navigation purposes.

This program is part of the WebTide software package.
/**********************************************************************/

/* PROGRAM: tidecor
   VERSION: 2.42
   DATE :   January 14, 2004
   Author : Jason Chaffey 
   Modifications: Shawn Oakey, Jason Chaffey
*/

/* Version 1.0 */
/*******************************************************************************
 * This program for tide correction of hydrographic data is derived from
 * interp2d.c - Modifications and additions made by Patrick Roussel May 21,1999
 ******************************************************************************/

/* Version 2.0 */
/* The FORTRAN components of Ver1.0 were translated to c for portability.
   Original FORTRAN written by M. Foreman at IOS.
   RAYBOUND bug fixed.
   Jason Chaffey              December 15, 2000
*/

/* Version 2.1 */
/* Added configuration file (tidecor2.1.cfg) so that the following are not
   hard-wired into the program:
       - mesh filenames (.nod and .ele);
       - the number of constituents to be used; and
       - for each constituent:
             - name of constituent; and
             - filename for data of that constituent.
   Jason Chaffey              January 15, 2001
*/

/* Version 2.2 */
/* Added support for complex input data (.v2c files)
   Jason Chaffey              March 23, 2001
*/

/* Version 2.21 */
/* Made various changes for use in WebTide
 *
 * Web(Tide/Drogue)
 *
 * Ocean Science Division
 * Bedford Institute of Oceanography
 * Fisheries and Oceans Canada
 *
 * Shawn Oakey
*/

/* Version 2.3 */
/*
 * Tested against Pawlawicz's T_Tide and 
 * Foreman's tide4_r2
 *
*/

/* Version 2.4 */
/*
 * Debugging (float - double and int - long int conflicts) done
 * Switched to Manhattan distance vs "real" distance in closestnode
 *          - speed increase consideration
 * Thanks to Herman Varma for the above
 * Jason Chaffey              August 7, 2003
 *
*/
/* Version 2.41 */
/*
 * Returned to using "real" distance calculations versus Manhattan distance.
 * Problems were discovered with Manhattan distance with the NE Pacific mesh.
 * Jason Chaffey              January 6, 2004
 *
*/
/* Version 2.42 */
/*
 * Bug in Manhattan distance calculation fixed. Returned to using it
 * for optimization.
 * Added fix for elements that cross International Dateline or
 * Greenwich Meridian in raybound.
 * Jason Chaffey              January 14, 2004
 *
*/
/* Version 2.5 */
/*
 * Fixed some problems with raybound for elements that cross International
 * Dateline or Greenwich Meridian in raybound. Problems became evident in the
 * global data set.
 * Jason Chaffey              June 16, 2008
 *
*/
/* Version 2.5.1 */
/*
 * Fixed numbering of nodes in RAYBOUND after old bug snuck back in.
 * Jason Chaffey              June 26, 2010
 *
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Define the value for 2 * pi */
#define twopi (2.0 * acos( -1.0 ))

/* This defines the (Julian) date of the switch to the Gregorian calendar */
/*   (Used in the julday subroutine) */
#define IGREG (15 + 31L * (10 + 12L * 1582 ))

/* These define the indices for the 5 constituents used in this program */
/* (Perhaps someday users will be allowed to choose which and how many to use */
#define m2  0
#define s2  1
#define n2  2
#define k1  3
#define o1  4

/* The following are structures used to hold the info of the constituents */
/* First is the structure to hold the satellites for each main constituent */
struct sattype
{ int deld[3];
  double phase;
  double ratio;
  int corr;
};
typedef struct sattype *satptr;

/* This is for the main constituents */
struct maintype
{ char *name;
  int dood[7];
  double phase;
  int nsats;
  satptr *sats;
  double f;
  double freq;
  double vu;
};
typedef struct maintype *mainptr;

/* This holds the main constituent factors for the shallow water constituents */
struct scontype
{ char *name;
  double factor;
};
typedef struct scontype *sconptr;

/* This structure holds the info for the shallow water constituents */
struct shalltype
{ char *name;
  int numcon;
  sconptr *shcon;
  double f;
  double freq;
  double vu;
};
typedef struct shalltype *shallptr;

/* Global constants: */
double R=6.3675e-8;         /* Radius of the earth */
long int nconst=5;          /* # of constituents */
int name_len=6;             /* Max. length of constituent names */
int maxcons = 45;           /* Initial max # of main constituents */
int maxshll = 110;          /* Initial max # of shallow water constituents */
    /* Program will re-alloc memory if more constituents loaded */

static int daytable[2][13] = {  /* Table of # days for each month */
  {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
  {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

/* Global variables: */
double *xdata, *ydata;	         /* Longitude and latitude */
double *xdata_base, *ydata_base; /* Base pointers for long and lat */
double *basis;                   /* Factors for interpolation */
int closest;                     /* Closest node number to an arbitrary point */
int *in, *in_base;               /* Pointers to the input data */
int numnodes, numeles;           /* Total # of nodes and elements */

/* Prototype declarations: */
satptr add_sat( char * );
void astro_angles( double, double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double * );
int basis2d( double, double);
char *caseup( char * );
int closestnode( int, double, double );
double dis_sq( double, double, double, double );
void get_date( int, int, long int *, long int * );
long int julday( long int, long int, long int );
int openvuf( char *, mainptr *, shallptr *, int *, int * );
double *phi2d( double *, double *, double, double );
int raybound( double *, double *, double, double );
void readneighbour( char * );
void readnodes( char *, char * );
void setvuf( long int, double, mainptr *, shallptr *, int, int );
double TideP( long int, long int, long int, long int, long int, double, double,
       double *, double *, char **, int, mainptr *, shallptr *, int, int );
int vuf( char *, mainptr *, shallptr *, int, int );

/*****************************************************************************/
satptr add_sat( char *satin )
/* Gets satellite info from input line, puts the info into a new
   satellite structure and returns the structure. */
{
  char *dummy;
  int d1, d2, d3;
  double phse, rat;
  satptr newsat;

  dummy = (char *)malloc( 3 * sizeof( char ));

  sscanf( satin, "%d %d %d %lf %lf", &d1,&d2,&d3, &phse, &rat );

  newsat = (satptr)malloc( sizeof( struct sattype ));

  newsat->deld[0] = d1;
  newsat->deld[1] = d2;
  newsat->deld[2] = d3;
  newsat->phase = phse;
  newsat->ratio = rat;

  if(( dummy = strchr( satin, 82 )) == NULL ) {
    newsat->corr = 0;
  } else {
    if( strncmp( dummy, "R1", 2 ) == 0 ) {
      newsat->corr = 1;
    } else if( strncmp( dummy, "R2", 2 ) == 0 ) {
      newsat->corr = 2;
    } else {
      newsat->corr = 0;
    }
  }
  return newsat;
}

/*****************************************************************************/
void astro_angles( double d1, double *h, double *pp, double *s, double *p,
     double *np, double *dh, double *dpp, double *ds, double *dp, double *dnp )
/* Calculates the following ephermides of the sun and moon:
   h  = mean longitude of the sun;
   pp = mean longitude of the solar perigee;
   s  = mean longitude of the moon;
   p  = mean longitude of the lunar perigee; and
   np = negative of the longitude of the mean ascending node.
   Also calculates their rate of change.
   Units are cycles ( cycles / 365 days for rates of change ).
   These formulae were taken from pp 98 and 107 of the Explanatory
   Supplement to the Astronomical Ephermeris and the American
   Ephermis and Nautical Almanac (1961) */
{
  double dum, d12, d2, d22, d23, f, f2;

  d12 = d1 * d1;
  d2 = d1 * 1.0E-04;
  d22 = d2 * d2;
  d23 = pow( d2, 3.0 );
  f = 360.0;
  f2 = f / 365.0;

  *h = ( 2.79696678E+02 + d1 * 9.856473354E-01 + d22 * 2.267E-05 ) / f;
  *h = modf( *h, &dum );

  *pp = ( 2.81220833E+02 + d1 * 4.70684E-05 + d22 * 3.39E-05 +
         d23 * 7.0E-08 ) / f;
  *pp = modf( *pp, &dum );

  *s = ( 2.70434164E+02 + d1 * 1.31763965268E+01 - d22 * 8.5E-05 +
         d23 * 3.9E-08 ) / f;
  *s = modf( *s, &dum );

  *p = ( 3.34329556E+02 + d1 * 1.114040803E-01 - d22 * 7.739E-04 -
         d23 * 2.6E-07 ) / f;
  *p = modf( *p, &dum );

  *np = ( -2.59183275E+02 + d1 * 5.29539222E-02 - d22 * 1.557E-04 -
         d23 * 5.0E-08 ) / f;
  *np = modf( *np, &dum );

  *dh = ( 9.856473354E-01 + d1 * 2.267E-05 * 2.0E-08 ) / f2;

  *dpp = ( 4.70684E-05 + d1 * 3.39E-05 * 2.0E-08 +
           d12 * 7.0E-08 * 3.0E-12 ) / f2;

  *ds = ( 1.31763965268E+01 - d1 * 8.5E-05 * 2.0E-08 +
           d12 * 3.9E-08 *3.0E-12 ) / f2;

  *dp = ( 1.114040803E-01 - d1 * 7.739E-04 * 2.0E-08 -
           d12 * 2.6E-07 * 3.0E-12 ) / f2;

  *dnp = ( 5.29539222E-02 - d1 * 1.557E-04 * 2.0E-08 -
           d12 * 5.0E-08 * 3.0E-12 ) / f2;
}

/*****************************************************************************/
int basis2d( double ptx, double pty )
/* Finds the closest node to a point (ptx, pty) and the element containing
   that point, if one exists.
   Also gets the basis functions for interpolations to that point.
   Returns the element number contaning the point (ptx, pty).
   Returns -999 if no element found containing the point (ptx, pty ). */
{
  int i, flag, ele, n11, n22, n33;
  double xlocal[3], ylocal[3];

  closest = closestnode( numnodes, ptx, pty );
  /*  Try to find an element that contains the point and the closest node. */
  ele = -999;
  in=in_base;
  for( i = 0; i <numeles; i++ ) {
    n11 = *in++;
    n22 = *in++;
    n33 = *in++;

    if(( closest == (n11) ) || ( closest == (n22) ) || ( closest == (n33) )) {
      xdata=xdata_base;
      xdata+=n11;
      ydata=ydata_base;
      ydata+=n11;
      xlocal[0] = *xdata;
      ylocal[0] = *ydata;

      xdata=xdata_base;
      xdata+=n22;
      ydata=ydata_base;
      ydata+=n22;
      xlocal[1] = *xdata;
      ylocal[1] = *ydata;

      xdata=xdata_base;
      xdata+=n33;
      ydata=ydata_base;
      ydata+=n33;
      xlocal[2] = *xdata;
      ylocal[2] = *ydata;
  /* See if the point is within this element */
      flag = raybound( &(xlocal[0]), &(ylocal[0]), ptx, pty );
      if( flag == 1 ) { /* The point is within the element */
        ele = i;
        basis = phi2d( &(xlocal[0]), &(ylocal[0]), ptx, pty );
        break;
      }
    }
  }

  /* if the closest node's elements don't work, search through all elements */
  if( ele < 0 ) {
    in=in_base;
    for( i = 0; i <numeles; i++ ) {
      xdata=xdata_base;
      xdata+=*in;
      xlocal[0] = *xdata;
      ydata=ydata_base;
      ydata+=*in;
      ylocal[0] = *ydata;
      in++;

      xdata=xdata_base;
      xdata+=*in;
      xlocal[1] = *xdata;
      ydata=ydata_base;
      ydata+=*in;
      ylocal[1] = *ydata;
      in++;

      xdata=xdata_base;
      xdata+=*in;
      xlocal[2] = *xdata;
      ydata=ydata_base;
      ydata+=*in;
      ylocal[2] = *ydata;
      in++;

  /* See if the point is within this element */
      flag = raybound( &(xlocal[0]), &(ylocal[0]), ptx, pty );
      if( flag == 1 ) {
        ele = i;
        basis = phi2d( &(xlocal[0]), &(ylocal[0]), ptx, pty );
        return( ele );
      }
    }
  }
  return( ele );
}

/*****************************************************************************/
char *caseup( char *lower )
/* Convert a string to all upper case. */
{
  char  toup[1], unique[50] = {0};
  int leng, i;
  
  leng = strlen( lower );
  for( i=0; i<leng; i++ ) {
    sscanf( lower, "%c", toup );;
    lower++;
	  unique[i] = toupper(toup[0]) ;
  }
  unique[leng] ='\0';

  return( strdup( unique ));
}

/*****************************************************************************/
int closestnode( int numnodes, double ptx, double pty )
/* Find the node that is closest to the point (ptx, pty). */
{
  int close, i;
  double currdist, closedist, d1,d2;

  xdata=xdata_base;
  ydata=ydata_base;
/*  closedist = dis_sq( *xdata, *ydata, ptx, pty );*/
/* Change to Manhattan distance - Aug 7 2003*/
  closedist = fabs( pty - *ydata ) + fabs( ptx - *xdata );
  xdata++;
  ydata++;
  close=0;

  for( i = 1; i < numnodes; i++ ) {
/*    currdist = dis_sq( *xdata, *ydata, ptx, pty );*/
/* Change to Manhattan distance - Aug 7 2003*/
    currdist = fabs( pty - *ydata ) + fabs( ptx - *xdata );
    xdata++;
    ydata++;
    if( currdist < closedist ) {
      closedist = currdist;
      close = i;
    }
  }
  return( close );
}

/*****************************************************************************/
double dis_sq( double lng1, double lat1, double lng2, double lat2 )
/* Calculate the distance between to latitude-longitude points */
{
  double dis_sqq;
  double deg_to_rad = asin(1.)/90.;

  dis_sqq = R * cos( 0.5 * deg_to_rad * ( lat1 + lat2 ));
  dis_sqq *= deg_to_rad * ( lng1 - lng2 );
  dis_sqq = dis_sqq * dis_sqq +
      R * R* deg_to_rad * deg_to_rad * ( lat1 - lat2 ) * ( lat1 - lat2 );
  return( dis_sqq );
}

/*****************************************************************************/
void get_date( int year, int dayofyear, long int *month, long int *day )
/* Get the day and month from the day # of the year */
{
  int i, leap;

  leap = (( year % 4 == 0 ) && ( year % 100 != 0 ) || ( year % 400 == 0 ));
  for( i = 1; dayofyear > daytable[leap][i]; i++ ) {
    dayofyear -= daytable[leap][i];
  }
  *month = (long)i;
  *day = (long)dayofyear;
}

/*****************************************************************************/
long int julday( long int id, long int im, long int iy )
/* Calculate the Julian day number.
   Accounts for the change to the Gregorian calandar. */
{
  long int jul, ja, jy, jm;

  jy = iy;
  if( jy == 0 ) {
    fprintf( stderr, "JULDAY: There is no year 0!\n" );
    exit( -1 );
  }
  if( jy < 0 ) ++jy;
  if( im > 2 ) {
    jm = im + 1;
  } else {
    --jy;
    jm = im + 13;
  }

  jul = (long int)( floor( 365.25 * jy ) + floor( 30.6001 * jm ) +
                id + 1720995 );
  if(( id + 31L * ( im + 12L * iy )) >= IGREG ) {
    ja = (long int)( 0.01 * jy );
    jul += 2 - ja + (long int)( 0.25 * ja );
  }
  return jul;
}

/*****************************************************************************/
/*****************************************************************************/
int main( int argc, char *argv[] )
{
  char *constfile, **constnames, *dummy, *elefile, *ext, *ext1, *nodefile;
  char *iosfile;
  FILE *cfg, **indatfile, *loctimefile, *tidecorfile;
  int cplx,i, j, ele, garbage, nmain, nread, nshall;
  long int day, month, year, hour, minute, dseconds, dayofyear, lgarbage;
  double *amp, *phase, *amp_base, *phase_base;
  double *amp2, *phase2, *amp_base2, *phase_base2;
  double longitude, latitude, reslt, reslt2, seconds;
  double elem_res[3],elem_res2[3], elem_min[3], elem_min2[3];
  double *min_tide, *min_tide_base;   /*Minimum tide values */
  double *min_tide2, *min_tide_base2; /*Minimum tide values (2nd output) */
  mainptr *cons;
  shallptr *shall;

  if(( argc < 3 ) || ( argc > 4)) {
    fprintf( stderr, "To start tide correction, " );
    fprintf( stderr, "enter the config, input and output files\n" );
    exit( 1 );
  }

  cplx = -9;

/* Read in parameters from the configuration file */
/* and open the tidal data files */

  if(( cfg = fopen( argv[1], "r" )) == NULL ) {
    fprintf( stderr, "Error in opening the cfg file!!\n" );
    exit( 9 );
  }

  dummy = (char *)malloc( 180 * sizeof( char ));
  nodefile = (char *)malloc( 180 * sizeof( char ));
  elefile =  (char *)malloc( 180 * sizeof( char ));
  iosfile =  (char *)malloc( 180 * sizeof( char ));

  fgets( nodefile, 180, cfg );
  strtok(nodefile,"\n");
  
  fgets(elefile, 180, cfg );
  strtok(elefile,"\n");

  /*Location of IOS Table*/
  fgets(iosfile, 180, cfg );
  strtok(iosfile,"\n");

  /* Read in the mesh */
  readnodes( nodefile, elefile );

  fgets( dummy, 180, cfg );
  j = sscanf( dummy, "%d", &nconst );
  if( j == EOF || j <= 0 ) {
    fprintf( stderr,"Error in reading the number of constituents.\n");
    fprintf( stderr, "Input line was %s\n", dummy );
    fprintf( stderr, "Cannot continue!\n\n" );
    exit( 9 );
  }

  constnames = (char **)malloc( nconst * sizeof( char * ));
  for( i = 0; i < nconst; i++ ) {
    constnames[i]=(char *)malloc( name_len * sizeof( char ));
  }
  indatfile = (FILE **)malloc( nconst * sizeof( FILE * ));
  constfile = ( char *)malloc( 180 * sizeof( char ));
  ext = (char *)malloc( 5 * sizeof( char ));
  ext1 = (char *)malloc( 5 * sizeof( char ));

  for( i = 0; i < nconst; i++ ) {
    fgets(constnames[i], 180, cfg );
    caseup( strtok(constnames[i],"\n") );
    fgets(constfile, 180, cfg );
    strtok(constfile,"\n");

    strcpy( ext, strrchr( constfile, 46 ) + 1 );
    ext = caseup( ext );
    if( cplx < 0 ) {
      if( strcmp( ext, "V2C" ) == 0 ) {
        cplx = 1;
      } else if( strcmp( ext, "S2C" ) == 0 ) {
        cplx = 0;
      } else if( strcmp( ext, "V2R" ) != 0 ) {
        fprintf( stderr, "Error. Unrecognized data file format (%s).\n", ext );
        fprintf( stderr, "Please check the filenames in the config file.\n" );
        exit( -12 );
      }
      strcpy( ext1, ext );
    } else {
      if( strcmp( ext, ext1 ) != 0 ) {
        fprintf(stderr, "ERROR! The extension for constituent data file %d", i);
        fprintf(stderr, " (%s)\n  is not the same as the first extension", ext);
        fprintf(stderr, " (%s).\n", ext1);
        fprintf( stderr, "Please check the filenames in the config file.\n" );
        exit( -13 );
      }
    }
    if(( cplx == 1 && strcmp( ext, "S2C" ) == 0 ) ||
       ( cplx == 0 && strcmp( ext, "V2C" ) == 0 )) {
      fprintf(stderr, "Error! Mixed types of constituent data files in the\n" );
      fprintf(stderr, "   configuration file.\nCannot Continue!\n" );
      exit( -14 );
    }
    if(( indatfile[i] = fopen( constfile, "r" )) == NULL ) {
      fprintf( stderr, "Error in opening constituent #%d datafile!!\n", i );
      fprintf( stderr, "Cannot continue!\n\n" );
     exit( 11 );
    } else {
      if( strcmp( constnames[i], "Z0" ) == 0 ) {
        fgets( dummy, 180, indatfile[i] );
        fgets( dummy, 180, indatfile[i] );
      } else {
        for( j = 0; j < 3; j++ ) {
          fgets( dummy, 180, indatfile[i] );
        }
      }
    }
  }
  if( i == nconst && cplx < 0 ) {
    fprintf(stderr, "Error.  Could not determine type of constituent data.\n");
    fprintf(stderr, " At least one extension must be either .v2c (currents)\n");
    fprintf(stderr, "   or .s2c (elevations).\n");
    fprintf(stderr, "Note: You cannot mix these types.\n" );
    fprintf(stderr, "Note: Z0 currents may be included with the .v2r" );
    fprintf(stderr, " extension.\n" );
    exit( -15 );
  }
  free( dummy );
  free( constfile );
  fclose( cfg );

/* Allocate memory */
  cons = (mainptr *)malloc( maxcons * sizeof( mainptr ));
  shall = (shallptr *)malloc( maxshll * sizeof( shallptr ));
  min_tide = (double *)malloc( numnodes * sizeof( double ));
  amp = (double *)malloc( numnodes * nconst * sizeof( double ));
  phase = (double *)malloc( numnodes * nconst * sizeof( double ));
  basis = (double *)malloc( 3 * sizeof( double ));
  if( cplx > 0 ) {
    amp2 = (double *)malloc( numnodes * nconst * sizeof( double ));
    phase2 = (double *)malloc( numnodes * nconst * sizeof( double ));
    min_tide2 = (double *)malloc( numnodes * sizeof( double ));
  }

/* Set pointer bases */
  min_tide_base = min_tide;
  amp_base = amp;
  phase_base = phase;
  if( cplx > 0 ) {
    min_tide_base2 = min_tide2;
    amp_base2 = amp2;
    phase_base2 = phase2;
  }

/* Load the model tidal data */
  for( j = 0; j < numnodes; j++ ) {
    for( i = 0; i < nconst; i++ ) {
      if( cplx == 0 ) {
        if(( nread = fscanf( indatfile[i], "%d %lf %lf",
                         &garbage, amp++, phase++)) != 3 ) {
          fprintf( stderr, "Error reading phase and amplitude file!!\n" );
          exit( 1 );
        }
      } else {
        if( strcmp( constnames[i], "Z0" ) == 0 ) {
          if(( nread = fscanf( indatfile[i], "%d %lf %lf",
                         &garbage, amp++, amp2++)) != 3 ) {
            fprintf(stderr, "Error reading mean current file %d!!\n", i );
            exit( 1 );
          }
          *phase = 0.0; phase++;
          *phase2 = 0.0; phase2++;
        } else {
          if(( nread = fscanf( indatfile[i], "%d %lf %lf %lf %lf",
                         &garbage, amp++, phase++, amp2++, phase2++ )) != 5 ) {
            fprintf( stderr,"%d %d %d\n", i, nread, garbage );
            fprintf(stderr, "Error reading vel phase and amplitude file %d!!\n",
                   i );
            exit( 1 );
          }
        }
      }
    }
  }

/* Calculate the minimum tides for each node */
  amp = amp_base;
  for( j = 0; j < numnodes; j++ ) {
    *min_tide = 0.0;
    if( cplx == 0 ) {
      for( i = 0; i < nconst; i++ ) {
/*        *min_tide = *min_tide + *amp++;*/
      }
    } else {
      *min_tide2 = 0.0;
      min_tide2++;
/*
      for( i = 0; i < nconst; i++ ) {
        *min_tide = *min_tide + *amp++;
        *min_tide2 = *min_tide2 + *amp2++;
      }
*/
    }
    min_tide++;
  }

/* Close the model tidal data files */
  for( i = 0; i < nconst; i++ ) {
    fclose( indatfile[i] );
  }
  free( indatfile );

/* Open the input/output files. */
  if(( loctimefile = fopen( argv[2], "r" )) == NULL ) {
    fprintf( stderr, "Error in opening the input file!!\n" );
    exit( 1 );
  }
  if(( tidecorfile = fopen( argv[3], "w" )) == NULL ) {
    fprintf( stderr, "Error in opening the output file!!\n" );
    exit( 1 );
  }

/* Set initial values. */
  day = 0;
  month = 0;
  year = 0;
  hour = 0;
  minute = 0;
  seconds = 0.0;
  latitude = 0.0;
  reslt = 0.0;
  reslt2 = 0.0;

/* Read in the constituent data */
  openvuf( iosfile, cons, shall, &nmain, &nshall );

/*  if(( cfg = fopen( "/home/chaffey/WebApp/test.txt", "w" )) == NULL ) {
    fprintf( stderr, "Error\n" );
    exit( 9 );
  }*/
/* Loop through the input file.
   For each line calculate the tidal correction */
  while(( nread = fscanf( loctimefile, "%lf %lf %ld %ld %ld %ld %lf",
                  &longitude, &latitude, &year, &dayofyear, &hour,
	          &minute, &seconds)) != EOF ) {
    if( nread != 7 ) {
      fprintf( stderr, "Error reading the input file!!\n" );
      exit( 1 );
    }
/*    fprintf( cfg, "IN: %lf %lf %ld %ld %ld %ld %lf\n",
      longitude, latitude, year, dayofyear, hour, minute, seconds );*/

    get_date( year, dayofyear, &month, &day );
/*    fprintf( stderr,"DATE: %d %d %d %d\n", year, month, day, dayofyear );*/
    ele = basis2d( longitude, latitude );

    if( ele < 0 ) {
/* Added for Web(Tide/Drogue) */
      fprintf(stderr,"Some Markers Not in the Domain\n");
      exit(-16);
/* Below is skipped for Web(Tide/Drogue) */
/* No element was found that contained the new position. */
/* Calculate the tidal correction for the closest node. */
      amp = amp_base;
      amp += closest*nconst;
      phase = phase_base;
      phase += closest*nconst;
      min_tide = min_tide_base;
      min_tide += closest;
      reslt = TideP( day, month, year, hour, minute, seconds, latitude,
             amp, phase, constnames, nconst, cons, shall, nmain, nshall );
      reslt = reslt + *min_tide;
      if( cplx > 0 ) {
        amp2 = amp_base2;
        amp2 += closest*nconst;
        phase2 = phase_base2;
        phase2 += closest*nconst;
        min_tide2 = min_tide_base2;
        min_tide2 += closest;
        reslt2 = TideP( day, month, year, hour, minute, seconds, latitude,
               amp2, phase2, constnames, nconst, cons, shall, nmain, nshall );
        reslt2 = reslt2 + *min_tide2;
      }

    } else {
/* If the point is inside and element, calculate the tidal correction for
   each node of the element and interpolate to get the tidal correction
   at the new position. */
      for( i = 0; i <= 2; i++ ) {
        in = in_base;
        in += i + ( 3 * ele );
        amp = amp_base;
        amp += (*in) * nconst;
        phase = phase_base;
        phase += (*in) * nconst;
        ydata = ydata_base;
        ydata += *in;
        xdata = xdata_base;
        xdata += *in;
        min_tide = min_tide_base;
        min_tide += *in;
        elem_res[i] = TideP( day, month, year, hour, minute, seconds, *ydata,
               amp, phase, constnames, nconst, cons, shall, nmain, nshall );
        elem_min[i] = *min_tide;
        if( cplx > 0 ) {
          amp2 = amp_base2;
          amp2 += (*in) * nconst;
          phase2 = phase_base2;
          phase2 += (*in) * nconst;
          min_tide2 = min_tide_base2;
          min_tide2 += *in;
          elem_res2[i] = TideP( day, month, year, hour, minute, seconds, *ydata,
               amp2, phase2, constnames, nconst, cons, shall, nmain, nshall );
          elem_min2[i] = *min_tide2;
        }
     }
     reslt = elem_res[0] * basis[0]
           + elem_res[1] * basis[1]
	   + elem_res[2] * basis[2]
	   + elem_min[0] * basis[0]
           + elem_min[1] * basis[1]
	   + elem_min[2] * basis[2];
     if( cplx > 0 ) {
       reslt2 = elem_res2[0] * basis[0]
              + elem_res2[1] * basis[1]
	      + elem_res2[2] * basis[2]
	      + elem_min2[0] * basis[0]
              + elem_min2[1] * basis[1]
	      + elem_min2[2] * basis[2];
      }
    }

/* Ouput the tidal correction */
    if( cplx == 0 ) {
    fprintf( tidecorfile,
          "%5.4lf %13.8lf %13.8lf %4ld %3ld %2ld %2ld %5.2lf\n", reslt,
          longitude, latitude, year, dayofyear, hour, minute, seconds);
    } else {
    fprintf( tidecorfile,
          "%5.4lf %5.4f %13.8lf %13.8lf %4ld %3ld %2ld %2ld %5.2lf\n", reslt,
          reslt2, longitude, latitude, year, dayofyear, hour, minute, seconds);
    }
  }

/*  fclose( cfg );*/
  fclose( loctimefile );
  fclose( tidecorfile );

  return( 0 );
}

/*****************************************************************************/
int openvuf( char *iosfile, mainptr *cons, shallptr *shall , int *nmain, int *nshall )
/* Read in the constituent data from the IOS_tidetbl file */
{
  FILE *VUF;
  int cnt, d1, d2, d3, d4, d5, d6, i, j, nln, nsat, num;
  double fact1, fact2, fact3, fact4, phse;
  char *dummy, *inp, *name1, *name2, *name3, *name4;
  char *satin1, *satin2, *satin3;
  mainptr newcon;
  shallptr newshall;
  sconptr newscon;

/* Allocate memory for input variables */
  inp = (char *)malloc( 92 * sizeof( char ));
  dummy = (char *)malloc( 21 * sizeof( char ));
  name1 = (char *)malloc( 21 * sizeof( char ));
  name2 = (char *)malloc( 21 * sizeof( char ));
  name3 = (char *)malloc( 21 * sizeof( char ));
  name4 = (char *)malloc( 21 * sizeof( char ));
  satin1 = (char *)malloc( 25 * sizeof( char ));
  satin2 = (char *)malloc( 25 * sizeof( char ));
  satin3 = (char *)malloc( 25 * sizeof( char ));

  if(( VUF = fopen( iosfile, "r" )) == NULL ) {
    fprintf( stderr, "Could not open >IOS_tidetbl<!!!\n" );
    fprintf( stderr, "Exiting ... \n\n" );
    exit( -11 );
  }

/* Counters for the number of main and shallow water constituents */
  *nmain = 0;
  *nshall = 0;

/* Read in the main constituents*/
  while( fgets( inp, 90, VUF ) != NULL ) {
    if( strlen( inp ) < (unsigned)name_len ) {
      break;
/* A blank line denotes the end of the main constituents */
    }

    sscanf( inp, "%s %d %d %d %d %d %d %lf %d", name1, &d1, &d2, &d3, &d4, &d5,
            &d6, &phse, &nsat );
    (*nmain)++;

/* Re-allocate memory if necessary */
    if( *nmain > maxcons ) {
      cons = (mainptr *)realloc( cons, *nmain * sizeof( mainptr ));
    }
    if( cons == NULL ) {
      fprintf(stderr," Error in allocating memory for main constituents.\n" );
      fprintf(stderr,  "Exiting ... \n\n" );
      exit( -21 );
    }

/* Create a new main constituent node */
    newcon = (mainptr)malloc( sizeof( struct maintype ));
    if( newcon ) {
      newcon->name = (char *)malloc( name_len * sizeof( char ));
      strncpy( newcon->name, name1, name_len );
      newcon->dood[0] = d1;
      newcon->dood[1] = d2;
      newcon->dood[2] = d3;
      newcon->dood[3] = d4;
      newcon->dood[4] = d5;
      newcon->dood[5] = d6;
      newcon->phase = phse;
      newcon->nsats = nsat;

/* Read in the satellites for this constituent, if any */
      if( nsat > 0 ) {
        j = 0;
        newcon->sats = (satptr *)malloc( nsat * sizeof( satptr ));
	nln = (( nsat - 1) / 3) + 1;
	for( i = 0; i < nln; i++ ) {
	  if( fgets( inp, 90, VUF ) == NULL ) {
	    fprintf(stderr,  "Error in reading %s satellite.\n", newcon->name );
	    fprintf(stderr,  "Exiting ... \n\n" );
	    exit( -12 );
	  }
	  cnt = nsat - ( i * 3 ); /* # of satellites on this line */
	  switch( cnt ) {
	    case 1:
	      sscanf( inp, "%12c%23c", dummy, satin1 );
	      newcon->sats[j] = add_sat( satin1 );
	      j++;
	      break;
	    case 2:
	      sscanf( inp, "%12c%23c%23c", dummy, satin1, satin2 );
	      newcon->sats[j] = add_sat( satin1 );
	      j++;
	      newcon->sats[j] = add_sat( satin2 );
	      j++;
	      break;
	    default:
	      sscanf( inp, "%12c%23c%23c%23c", dummy, satin1, satin2, satin3 );
	      newcon->sats[j] = add_sat( satin1 );
	      j++;
	      newcon->sats[j] = add_sat( satin2 );
	      j++;
	      newcon->sats[j] = add_sat( satin3 );
	      j++;
	      break;
	  }
	}
      } else {
        /*  NO satellites */
        newcon->sats = NULL;
      }
      cons[*nmain-1] = newcon;
    }
  }

/* Read in the shallow water constiuents */
  while( fgets( inp, 90, VUF ) != NULL ) {
    if( strlen( inp ) < (unsigned)name_len ) {
/* A blank line denotes the end of the shallow water constituents */
      break;
    }

    sscanf( inp, "%s %d", name1, &num );
    (*nshall)++;
/* Re-allocate memory if necessary */
    if( *nshall > maxshll ) {
      shall = (shallptr *)realloc( shall, *nshall * sizeof( shallptr ));
    }
    if( shall == NULL ) {
      fprintf(stderr," Error in allocating memory for shallow water const.\n" );
      fprintf(stderr, "Exiting ... \n\n" );
      exit( -21 );
    }

/* Create a new shallow water constituent node */
    newshall = (shallptr)malloc( sizeof( struct shalltype ));
    if( newshall ) {
      newshall->name = (char *)malloc( name_len * sizeof( char ));
      strncpy( newshall->name, name1, name_len );
      strcpy( name1, "" ); strcpy( name2, "" );
      strcpy( name3, "" ); strcpy( name4, "" );
      newshall->numcon = num;
      newshall->shcon = (sconptr *)malloc( num * sizeof( sconptr ));
      switch( num ) {  /* # of main constituent factors */
/* For each factor, create a node, get the info and add it to the shallow
   water constituent node */
        case 4:
          sscanf( inp, "%s %d %lf %s %lf %s %lf %s %lf %s", dummy, &i, &fact1,
	        name1, &fact2, name2, &fact3, name3, &fact4, name4 );
	  newscon = (sconptr)malloc( sizeof( struct scontype ));
          newscon->name = (char *)malloc( name_len * sizeof( char ));
	  strncpy( newscon->name, name4, name_len );
	  newscon->factor = fact4;
	  newshall->shcon[3] = newscon;
        case 3:
          sscanf( inp, "%s %d %lf %s %lf %s %lf %s", dummy, &i, &fact1, name1,
	        &fact2, name2, &fact3, name3 );
	  newscon = (sconptr)malloc( sizeof( struct scontype ));
          newscon->name = (char *)malloc( name_len * sizeof( char ));
	  strncpy( newscon->name, name3, name_len );
	  newscon->factor = fact3;
	  newshall->shcon[2] = newscon;
        case 2:
          if( strlen( name2 ) == 0 ) {
            sscanf( inp, "%s %d %lf %s %lf %s", dummy, &i, &fact1, name1,
	        &fact2, name2 );
	  }
	  newscon = (sconptr)malloc( sizeof( struct scontype ));
          newscon->name = (char *)malloc( name_len * sizeof( char ));
	  strncpy( newscon->name, name2, name_len );
	  newscon->factor = fact2;
	  newshall->shcon[1] = newscon;
        case 1:
          if( strlen( name1 ) == 0 ) {
            sscanf( inp, "%s %d %lf %s", dummy, &i, &fact1, name1 );
	  }
	  newscon = (sconptr)malloc( sizeof( struct scontype ));
          newscon->name = (char *)malloc( name_len * sizeof( char ));
	  strncpy( newscon->name, name1, name_len );
	  newscon->factor = fact1;
	  newshall->shcon[0] = newscon;
          break;
       }
      shall[*nshall-1] = newshall;
    }
  }

  fclose( VUF );

/* Free the memory of the input variables */
  free( inp );
  free( dummy );
  free( name1 );
  free( satin1 );
  free( name2 );
  free( satin2 );
  free( name3 );
  free( satin3 );
  free( name4 );

  return( 0 );
}

/*****************************************************************************/
double *phi2d( double *xloc, double *yloc, double ptx, double pty )
/* Calculates the basis functions for interpolating to a point inside
   an element. */
{
  double *phi, area, a, b, c;
  int i, j, k;

  phi = (double *)malloc( 3 * sizeof( double ));

  area = 0.5 * ( xloc[0] * ( yloc[1] - yloc[2] ) +
                xloc[1] * ( yloc[2] - yloc[0] ) +
                xloc[2] * ( yloc[0] - yloc[1] ));

  /* Calculate the Basis function... */
  for( i = 0; i <= 2; i++ ) {
    switch( i ) {
      case 0 : j = 1;  k = 2;  break;
      case 1 : j = 2;  k = 0;  break;
      case 2 : j = 0;  k = 1;  break;
    }
    a = ( xloc[j] * yloc[k] - xloc[k] * yloc[j] ) / ( area * 2 );
    b = ( yloc[j] - yloc[k] ) / ( area * 2 );
    c = -1 * ( xloc[j] - xloc[k] ) / ( area * 2 );
    phi[i] = a + b * ptx + c * pty;
  }
  return( phi );
}


/******************************************************************************/
int raybound( double *xd, double *yd, double ptx, double pty )
/*  Subroutine to check wether or not a point is inside a polygon.
The process is as follows:
	Use an arbitrary ray (here, y = constant and x >= xref), starting from 
the point and going off to infinity.
	Count the number of polygon boundaries it crosses.
	If an odd number, the point is inside the polygon, otherwise it is
outside.   */
{
  int i, j, bcross;
  double b, m, x;
  
  bcross = 0; /* Number of boundary crossings. */

/* Check to see if the element side crosses the International Dateline
   (changes sign at +180/-180 degrees) and if so, change the longitudes
   so that they all have the same sign. */
  
  if( ptx > 0.0 ) { 
    if(( xd[0] < -170.0 ) && ((xd[1] > 170.0 ) || ( xd[2] > 170.0 ))) {
      xd[0] += 360.0;
    }
    if(( xd[1] < -170.0 ) && ((xd[0] > 170.0 ) || ( xd[2] > 170.0 ))) {
      xd[1] += 360.0;
    }
    if(( xd[2] < -170.0 ) && ((xd[1] > 170.0 ) || ( xd[0] > 170.0 ))) {
      xd[2] += 360.0;
    }
  } else {
    if(( xd[0] > 170.0 ) && ((xd[1] < -170.0 ) || ( xd[2] < -170.0 ))) {
      xd[0] -= 360.0;
    }
    if(( xd[1] > 170.0 ) && ((xd[0] < -170.0 ) || ( xd[2] < -170.0 ))) {
      xd[1] -= 360.0;
    }
    if(( xd[2] > 170.0 ) && ((xd[1] < -170.0 ) || ( xd[0] < -170.0 ))) {
      xd[2] -= 360.0;
    }
  }
  
/* As above, except for the Greenwich meridian, for longitude coordinates
   that range from 0 to 360 degrees. */

  if( ptx > 350.0 ) {
    if(( xd[0] < 10.0 ) && ((xd[1] > 350.0 ) || ( xd[2] > 350.0 ))) {
      xd[0] += 360.0;
    }
    if(( xd[1] < 10.0 ) && ((xd[0] > 350.0 ) || ( xd[2] > 350.0 ))) {
      xd[1] += 360.0;
    }
    if(( xd[2] < 10.0 ) && ((xd[1] > 350.0 ) || ( xd[0] > 350.0 ))) {
      xd[2] += 360.0;
    }
  } else {
    if(( xd[0] > 350.0 ) && ((xd[1] < 10.0 ) || ( xd[2] < 10.0 ))) {
      xd[0] -= 360.0;
    }
    if(( xd[1] > 350.0 ) && ((xd[0] < 10.0 ) || ( xd[2] < 10.0 ))) {
      xd[1] -= 360.0;
    }
    if(( xd[2] > 350.0 ) && ((xd[1] < 10.0 ) || ( xd[0] < 10.0 ))) {
      xd[2] -= 360.0;
    }
  }
  
  for( i = 0; i <= 2; i++ ) {
  
  /* for each line segment around the element */
    j = (( i == 2 ) ? 0 : i + 1 );
   
  /* If both endpoints of the line segment are on the same (vertical)
side of the ray, do nothing.
     Otherwise, count the number of times the ray intersects the segment. */
    if( !((( yd[i] < pty ) && ( yd[j] < pty )) ||
         (( yd[i] >= pty ) && ( yd[j] >= pty )))) {

      if( xd[i] != xd[j] ) {
        m = ( yd[j] - yd[i] ) / ( xd[j] - xd[i] );
        b = yd[i] - m * xd[i] ;
        x = ( pty - b ) / m ;
        if( x > ptx ) { 
          bcross++;
        }
      } else {
        if( xd[j] > ptx ) {
          bcross++;
        }
      }
    }
  }

/*  Return the evenness/oddness of the boundary crossings
        i.e. the remainder from division by two. */
  return( bcross % 2 );
}

/*****************************************************************************/
void readnodes( char *nfile, char *efile )
/* Read in the mesh from a .nod/.ele set of files. */
{
  FILE *nodal, *eltal, *WC;
  double rdum, sdum;
  int i, in1, in2, in3, dum, nnd, nne, zero;
  char wccom[79], rmcom[79];

  /* Read Node file */
  nodal = fopen( nfile, "r" );
  fscanf( nodal, "%d %lf %lf", &dum, &rdum, &sdum );
/*  fprintf( stdout, "%d %lf %lf", dum, rdum, sdum );*/
  numnodes = 0;
  while( !feof( nodal )) {
    numnodes++;
    fscanf( nodal, "%d %lf %lf", &dum, &rdum, &sdum );
  }

  /* Read Element file */
  eltal = fopen( efile, "r" );
  fscanf( eltal, "%d %d %d %d", &dum, &in1, &in2, &in3 );
  numeles = 0;
  while( !feof( eltal )) {
    numeles++;
    fscanf( eltal, "%d %d %d %d", &dum, &in1, &in2, &in3 );
  }

  /* Read in the nodes. */
  rewind( nodal );
  xdata = (double *)malloc( (numnodes) * sizeof( double ));
  ydata = (double *)malloc( (numnodes) * sizeof( double ));
  xdata_base = xdata;
  ydata_base = ydata;
  for( i = 0; i < numnodes; i++ ) {
    fscanf( nodal, "%d %lf %lf", &dum, xdata++, ydata++ );
/* The program starts counting at zero. Does the data file? */
    if( i == 0 ) {
      zero = ( dum == 0 ) ? 0 : 1;
    }
 }
  fclose( nodal );

  /* Read in the elements. */
  rewind( eltal );
  in = (int *)malloc(3 * numeles * sizeof( int ));
  in_base = in;
  for( i = 0; i <numeles; i++ ) {
    fscanf( eltal, "%d %d %d %d", &dum, &in1, &in2, &in3 );
    *in++ = in1;
    *in++ = in2;
    *in++ = in3;
  }
  if( zero != 0 ) {
/* if the node file started counting nodes at 1 instead of zero,
   must decrement node #'s in the element list to correspond to the
   program counting from zero. */
    in=in_base;
    for( i = 0; i < numeles; i++ ) {
      in[i*3]--;
      in[i*3+1]--;
      in[i*3+2]--;
    }
  }
  fclose( eltal );
}

/*****************************************************************************/
void setvuf(long int kh, double xlat, mainptr *cons, shallptr *shall, int nmain,
             int nshall )
/* Calculate the amplitudes, phases, etc. for each of the constituents */
{
  double d1, hh, tau, dtau, slat, sumc, sums, v, vdbl;
  double adj, dum, dd[6], uu, uudbl;
  double h, pp, s, p, enp, dh, dpp, ds, dp, dnp;
  int i, j, k,indx;
  long int kd, ktmp;
  satptr sat;
/*  FILE *jtest;*/

  kd = julday( 31, 12, 1899 );
  d1 = (double)(kh - kd) - 0.5;
  ktmp=kh*24;

  astro_angles( d1, &h, &pp, &s, &p, &enp, &dh, &dpp, &ds, &dp, &dnp );
  hh = (double)ktmp - ( floor( (double)( ktmp / 24.0 )) * 24.0 );
  tau = hh / 24.0 + h - s;
  dtau = 365.00 + dh - ds;
  slat = sin(( twopi / 2.0 ) * ( xlat / 180.0 ));

/* The main constituents */
  for( k = 0; k < nmain; k++ ) {
    for( i = 0; i < 6; i++ ) {
      dd[i] = cons[k]->dood[i];
    }
    cons[k]->freq = ( dd[0] * dtau + dd[1] * ds + dd[2] * dh + dd[3] *dp +
                      dd[4] * dnp + dd[5] * dpp ) / ( 24.0 * 365.0 );
    vdbl = dd[0] * tau + dd[1] * s + dd[2] * h + dd[3] * p +
                      dd[4] * enp + dd[5] * pp + cons[k]->phase;
    v = vdbl - ( floor( floor( vdbl ) / 2.0 ) * 2.0 );

    sumc = 1.0;
    sums = 0.0;

    for( i = 0; i < cons[k]->nsats; i++ ) {
      sat = cons[k]->sats[i];
      switch( sat->corr ) {
        case 0:
	  adj = sat->ratio;
	  break;
	case 1:
	  adj = sat->ratio * 0.36309 * ( 1.0 - 5.0 * slat * slat ) / slat;
	  break;
	case 2:
	  adj = sat->ratio * 2.59808 * slat;
	  break;
      }
      uudbl = (double)sat->deld[0] * p + (double)sat->deld[1] * enp +
              (double)sat->deld[2] * pp + (double)sat->phase;
      uu = modf( uudbl, &dum );
      sumc += ( adj * cos( uu * twopi ));
      sums += ( adj * sin( uu * twopi ));
    }
    cons[k]->f = sqrt(( sumc * sumc ) + ( sums * sums ));
    cons[k]->vu = v + atan2( sums, sumc ) / twopi;
  }

/* The shallow water constituents */
  for( k = 0; k < nshall; k++ ) {
    shall[k]->f = 1.0;
    shall[k]->vu = 0.0;
    shall[k]->freq = 0.0;
    for( i = 0; i < shall[k]->numcon; i++ ) {
      for( j = 0; j < nmain; j++ ) {
        if( strcmp( cons[j]->name, shall[k]->shcon[i]->name ) == 0 ) {
          shall[k]->f *= pow( cons[j]->f, fabs( shall[k]->shcon[i]->factor ));
	  shall[k]->vu += ( shall[k]->shcon[i]->factor * cons[j]->vu );
	  shall[k]->freq += ( shall[k]->shcon[i]->factor * cons[j]->freq );
	  break;
	}
      }
    }
  }
/*  indx = vuf( "K1", cons, shall, nmain, nshall );
  if((jtest = fopen( "/home/chaffey/jtest.dat", "a+" )) == NULL ) {
    fprintf( stderr, "Error opening test file.\n" );
    exit(109);
  }
  fprintf( jtest, "setvuf : %f %f\n", cons[indx]->f,cons[indx]->vu );
  fclose( jtest );*/
}

/*****************************************************************************/
double TideP( long int day, long int month, long int year, long int hour, 
       long int minute, double seconds, double latitude, double *ampl,
       double *phase, char **innames, int numnames, mainptr *cons,
       shallptr *shall, int nmain, int nshall )
/* Calculates and returns the tidal correction */
{
/*  FILE *jtest;*/
  int i, indx;
  long int kd;
  double dthr, dum, radgmt, revgmt, res;

  kd = julday( day, month, year );

  setvuf( kd, latitude, cons, shall, nmain, nshall );
  dthr = (( (double)hour * 3600.0 ) + ( (double)minute * 60.0 ) +
            (double)seconds ) / 3600.0;
  res = 0.0;

/*  if((jtest = fopen( "/home/chaffey/jtest.dat", "a+" )) == NULL ) {
    fprintf( stderr, "Error opening test file.\n" );
    exit(109);
  }*/
/*  fprintf( jtest, "%d %d %d %d %d %f\n",day,month,year,hour,minute,seconds);*/
/* For each of the desired constituents ... (See top of program) */
  for( i = 0; i < numnames; i ++ ) {
/* Find the constituent from those loaded from IOS_tidetbl */
    indx = vuf( innames[i], cons, shall, nmain, nshall );
    if( indx < 0 ) {
      fprintf( stderr, "Bad Input Constituent: %s\n", innames[i] );
      fprintf( stderr, "Exiting ... \n\n" );
      exit( -1 );
    }

    if( indx < nmain ) {                        /* Main constituent */
      revgmt = cons[indx]->freq * dthr + cons[indx]->vu - phase[i] / 360.0;
      radgmt = twopi * modf( revgmt, &dum );
      res += cons[indx]->f * ampl[i] * cos( radgmt );
/*      fprintf( jtest, "%f %f %f %f %f %f\n", cons[indx]->freq, dthr, cons[indx]->vu,
         phase[i],radgmt, cons[indx]->f );*/
    } else if(( indx - nmain ) < nshall ) {     /* Shallow water constituent */
      indx -= nmain;
      revgmt = shall[indx]->freq * dthr + shall[indx]->vu - phase[i] / 360.0;
      radgmt = twopi * modf( revgmt, &dum );
      res += shall[indx]->f * ampl[i] * cos( radgmt );

    } else {
      fprintf( stderr, "Error in index: %d %d %d\n", indx, nmain, nshall );
      fprintf( stderr, "Exiting ... \n\n" );
      exit( -2 );
    }
  }
/*  fclose( jtest );*/
  return( res );
}

/*****************************************************************************/
int vuf( char *inname, mainptr *cons, shallptr *shall, int nmain, int nshall )
/* Finds constituent info corresponding to inname and returns the index to
   the node containing the info. Shallow water constituent indices are
   returned as their number greater than the max # of main constituents.
   e.g. if we want the 2nd shallow water constituent and there are 45
   main constituents, then the indice returned is 46, since constituents
   are counted from zero. ( 45 - 1 + 2 = 46 ) */
{
  int i, j;

  i = 0;
  while( strcmp( cons[i]->name, inname ) != 0 ) {
    i++;
    if( i == nmain ) break;
  }
  if( i == nmain ) {
    j = 0;
    while( strcmp( shall[j]->name, inname ) != 0 ) {
      j++;
      if( j == nshall ) break;
    }
  }

  if( i < nmain ) {
    return( i );
  } else if( j < nshall ) {
    return( i + j );
  } else {
    fprintf( stderr, "Constituent %s not found!\n", inname );
    fprintf( stderr, "Exiting ... \n\n" );
    exit( -9 );
  }
}
