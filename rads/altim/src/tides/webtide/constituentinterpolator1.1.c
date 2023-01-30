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

/* PROGRAM: constituentinterpolator
   VERSION: 1.0
   DATE :   Nov 10, 2007
   Author : Chad Gilbert 
   Modifications: Shawn Oakey, Jason Chaffey
*/

/* Version 1.0 */
/*******************************************************************************
 * This program for use with WebTide is derived from tidecor
 ******************************************************************************/
/* Version 1.1 */
/*******************************************************************************
 * Fixed some problems with raybound for elements that cross International
 * Dateline or Greenwich Meridian in raybound. Problems became evident in the
 * global data set.
 * Jason Chaffey              June 26, 2010
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Define the value for 2 * pi */
#define twopi (2.0 * acos( -1.0 ))

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
int basis2d( double, double);
char *caseup( char * );
int closestnode( int, double, double );
double *phi2d( double *, double *, double, double );
int raybound( double *, double *, double, double );
void readnodes( char *, char * );

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
/*****************************************************************************/
int main( int argc, char *argv[] )
{
  char *constfile, **constnames, *dummy, *elefile, *ext, *ext1, *nodefile;
  char *iosfile;
  FILE *cfg, **indatfile, *loctimefile, *tidecorfile;
  int cplx, i, j, ele, garbage, nmain, nread;
  long int day, month, year, hour, minute, dseconds, dayofyear, lgarbage;
  double *amp, *phase, *amp_base, *phase_base;
  double *amp2, *phase2, *amp_base2, *phase_base2;
  double longitude, latitude, reslt, reslt2, seconds,latold = 0, longold = 0;
  mainptr *cons;
/* added Jan 8 2007 for Chad's Get Constituent application */
  double **Camp, **Camp2, **Cphase, **Cphase2;
  double *Cresult_amp, *Cresult_amp2, *Cresult_phase, *Cresult_phase2;
  double *X, *Y, *X2, *Y2;
/* end added Jan 8 2007 */

  if(( argc < 3 ) || ( argc > 4)) {
    fprintf( stderr, "To start tide correction, " );
    fprintf( stderr, "enter the config, input and output files\n" );
    exit( 1 );
  }

  cplx = -9;

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

  fgets(iosfile, 180, cfg );
  strtok(iosfile,"\n");

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

  cons = (mainptr *)malloc( maxcons * sizeof( mainptr ));
  amp = (double *)malloc( numnodes * nconst * sizeof( double ));
  phase = (double *)malloc( numnodes * nconst * sizeof( double ));
  basis = (double *)malloc( 3 * sizeof( double ));
  if( cplx > 0 ) {
    amp2 = (double *)malloc( numnodes * nconst * sizeof( double ));
    phase2 = (double *)malloc( numnodes * nconst * sizeof( double ));
  }

  amp_base = amp;
  phase_base = phase;
  if( cplx > 0 ) {
    amp_base2 = amp2;
    phase_base2 = phase2;
  }

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
            fprintf(stderr, "Error reading vel phase and amplitude file %d!!\n", i );
            exit( 1 );
          }
        }
      }
    }
  }

  amp = amp_base;

  for( i = 0; i < nconst; i++ ) {
    fclose( indatfile[i] );
  }
  free( indatfile );

  if(( loctimefile = fopen( argv[2], "r" )) == NULL ) {
    fprintf( stderr, "Error in opening the input file!!\n" );
    exit( 1 );
  }
  if(( tidecorfile = fopen( argv[3], "w" )) == NULL ) {
    fprintf( stderr, "Error in opening the output file!!\n" );
    exit( 1 );
  }

  Camp = (double **)malloc( nconst * sizeof( double * ));
  Camp2 = (double **)malloc( nconst * sizeof( double * ));
  Cphase = (double **)malloc( nconst * sizeof( double * ));
  Cphase2 = (double **)malloc( nconst * sizeof( double * ));
  for( i = 0; i < nconst; i++ ) {
    Camp[i] = (double *)malloc( 3 * sizeof( double ));
    Camp2[i] = (double *)malloc( 3 * sizeof( double ));
    Cphase[i] = (double *)malloc( 3 * sizeof( double ));
    Cphase2[i] = (double *)malloc( 3 * sizeof( double ));
  }
  Cresult_amp = (double *)malloc( nconst * sizeof( double ));
  Cresult_amp2 = (double *)malloc( nconst * sizeof( double ));
  Cresult_phase = (double *)malloc( nconst * sizeof( double ));
  Cresult_phase2 = (double *)malloc( nconst * sizeof( double ));
  X = (double *)malloc( nconst * sizeof( double ));
  Y = (double *)malloc( nconst * sizeof( double ));
  X2 = (double *)malloc( nconst * sizeof( double ));
  Y2 = (double *)malloc( nconst * sizeof( double ));
  
/*  double Camp[nconst][3];
  double Camp2[nconst][3];
  double Cphase[nconst][3];
  double Cphase2[nconst][3];
  double Cresult_amp[nconst];
  double Cresult_amp2[nconst];
  double Cresult_phase[nconst];
  double Cresult_phase2[nconst];
  double X[nconst];
  double Y[nconst];
  double X2[nconst];
  double Y2[nconst];
  double latold = 0;
  double longold = 0;*/

  while(( nread = fscanf( loctimefile, "%lf %lf %ld %ld %ld %ld %lf",
                  &longitude, &latitude, &year, &dayofyear, &hour,
	          &minute, &seconds)) != EOF ) {

    if(longitude == longold && latitude == latold) {
      continue;
    }
    longold = longitude;
    latold = latitude;
    
/*    fprintf( stderr, "READ: %f %f\n", longitude, latitude );*/
	
    ele = basis2d(longitude, latitude);

/*    fprintf( stderr, "ele is %d\n", ele );*/

    for( i = 0 ; i < nconst ; i++ ) {
      for( j = 0 ; j < 3 ; j++ ) {
        in = in_base + j + ele*3;
        amp = amp_base + i + (*in) * nconst;
        phase = phase_base + i + (*in) * nconst;
        Camp[i][j] = *amp;
        Cphase[i][j] = *phase*twopi/360.0;
      }
    }

/*    fprintf( stderr, "Before cplx.\n" );*/
    
    if( cplx > 0 ) {
      for( i = 0 ; i < nconst ; i++ ) {
        for( j = 0 ; j < 3 ; j++ ) {
          in = in_base + j + ele*3;
          amp2 = amp_base2 + i + (*in) * nconst;
          phase2 = phase_base2 + i + (*in) * nconst;
          Camp2[i][j] = *amp2;
          Cphase2[i][j] = *phase2*twopi/360.0;
        }
      }
    }
/*    fprintf( stderr, "After cplx.\n" );*/

    for( i = 0 ; i < nconst ; i++ ) {
      X[i] = Camp[i][0] * cos(Cphase[i][0]) * basis[0] +
             Camp[i][1] * cos(Cphase[i][1]) * basis[1] +
             Camp[i][2] * cos(Cphase[i][2]) * basis[2];
		
      Y[i] = Camp[i][0] * sin(Cphase[i][0]) * basis[0] +
             Camp[i][1] * sin(Cphase[i][1]) * basis[1] +
             Camp[i][2] * sin(Cphase[i][2]) * basis[2];

      Cresult_amp[i] = sqrt(X[i]*X[i] + Y[i]*Y[i]);
      Cresult_phase[i] = atan2(Y[i],X[i]) * 360.0/twopi;
      if(Cresult_phase[i] < 0) {
        Cresult_phase[i] += 360.0;
      }
    }

/*    fprintf( stderr, "Scalar results.\n" );*/
    
    if( cplx > 0 ) {
      for( i = 0 ; i < nconst ; i++ ) {
        X2[i] = Camp2[i][0] * cos(Cphase2[i][0]) * basis[0] +
                Camp2[i][1] * cos(Cphase2[i][1]) * basis[1] +
               Camp2[i][2] * cos(Cphase2[i][2]) * basis[2];
		  
        Y2[i] = Camp2[i][0] * sin(Cphase2[i][0]) * basis[0] +
                Camp2[i][1] * sin(Cphase2[i][1]) * basis[1] +
                Camp2[i][2] * sin(Cphase2[i][2]) * basis[2];

        Cresult_amp2[i] = sqrt(X2[i]*X2[i] + Y2[i]*Y2[i]);
        Cresult_phase2[i] = atan2(Y2[i],X2[i]) * 360.0/twopi;
        if( Cresult_phase2[i] < 0 ) {
          Cresult_phase2[i] += 360.0;
        }
      }
    }
/*    fprintf( stderr, "CPLX results.\n" );*/

    if( cplx == 0 ) {
      for( i = 0 ; i < nconst ; i++ ) {
        fprintf( tidecorfile, "%s:\t%lf\t%lf\t%lf\t%3.6lf\n", constnames[i],
                 longitude, latitude, Cresult_amp[i], Cresult_phase[i] );
      }
    } else {
      for( i = 0 ; i < nconst ; i++ ) {
        fprintf( tidecorfile, "%s:\t%lf\t%lf\t%lf\t%lf\t%lf\t%3.6lf\n",
                 constnames[i], longitude, latitude, Cresult_amp[i],
                 Cresult_phase[i], Cresult_amp2[i], Cresult_phase2[i] );
      }
    }
    fprintf( tidecorfile, "\n");
  }


  fclose( loctimefile );
  fclose( tidecorfile );

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
        if( x > ptx ) bcross++;
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
