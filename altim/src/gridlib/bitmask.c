/*BITMASK -- Routines to load and query geographical bit masks
 +
      SUBROUTINE BITMASKINIT (BITMASK)
      CHARACTER*(*) BITMASK

      FUNCTION BITMASKEVAL (LAT, LON)
      REAL*8  LAT, LON
      INTEGER BITMASKEVAL

      SUBROUTINE BITMASKEND ()

  The routine BITMASKINIT reads the bitmap BITMASK into memory.
  This step is manditory before BITMASKEVAL can be used. BITMASK is
  a binary PBM (portable bit map) file that can be viewed with programs
  like XV. Generally we choose a bit status 0 (black) for land and
  1 (white) for ocean. The bitmap has to cover the area 0 to 360 degrees
  longitude and -90 to 90 degrees latitude and should have an equal
  number of bits per degree in both directions. Therefore, the mask
  should be twice as wide as it is high.

  The integer function BITMASKEVAL returns 1 if the bitmap indicates
  land at the location (LAT, LON), 0 if the location is ocean, sea, or
  lake (water), -1 if an error occurred. LAT is latitude between -90 and 90
  degrees, LON may be given in an "unlimited" range.

  The routine BITMASKEND removes the BITMASK from memory and frees
  it for further usage.

  Arguments:
    BITMASK      (input): Name of the bitmap (PBM) file that has the mask
    LAT          (input): Latitude (degrees)
    LON          (input): Longitude (degrees) (range unlimited)
    BITMASKEVAL (output): 1=land, 0=water, -1=error
--
     Feb-1999 - Created by Ejo Schrama
  03-Mar-1999 - Adapted by Remko Scharroo
------------------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "machine.h"
#include <sysdep.h>
#ifdef CAPITALS
#define bitmaskinit BITMASKINIT
#define bitmaskeval BITMASKEVAL
#define bitmaskend BITMASKEND
#endif
#ifdef UNDERSCOREAFTER
#define bitmaskinit bitmaskinit_
#define bitmaskeval bitmaskeval_
#define bitmaskend bitmaskend_
#endif
#ifdef UNDERSCOREBEFORE
#define bitmaskinit _bitmaskinit
#define bitmaskeval _bitmaskeval
#define bitmaskend _bitmaskend
#endif
//
// some global variables in the problem
// ------------------------------------
//
unsigned char **buffer;
unsigned char *maskset;
int4 nrow,ncol;
real8 bpd; 
//
// initialize the bitmask
// ----------------------
// 
void bitmaskinit( char *inputfile, int slen )
{
  int i,nstr,nret,pos;
  int ncol8;
  int n;
  FILE *fp;
  char str[160];
  // Copy the F77 filename to temporary string
  i = slen - 1;
  strncpy(str,inputfile,slen);
  while (str[i] == ' ' && i >= 0) i--;
  str[i+1]='\0';
  //
  // Open the file
  // -------------
  //
  //printf("Opening file %s\n",str);
  fp = fopen(str,"r");
  if (fp == NULL) { printf("Can't find file %s\n",str); exit(-1); } 
  //
  // Print the three header lines
  // ----------------------------
  //
  nstr = 0;
  pos = 0;
  nret = 0;
  //printf("***** File header\n\n");
  while ((nret < 3) && (!feof(fp))) {
    pos++;
    if ((str[nstr++] = fgetc(fp)) == 10) {
      str[nstr] = 0;
      //printf("%s",str);
      nret++;
      nstr=0;
    }
  }
  //printf("\n***** Currently at pos %d\n",pos);
  //
  // # rows and columns in this file (it is on the last record produced by gs)
  // -------------------------------------------------------------------------
  //
  sscanf(str,"%d %d",&ncol,&nrow);
  if (ncol != 2*nrow) {
    printf("Your mask file has the wrong dimensions!\n");
    exit(-1);
  }
  //
  // since 8 bits are stored in a character, ncol / 8 is the amount of char's
  // this will automatically take care of padding bits at the end of the record
  // --------------------------------------------------------------------------
  //
  ncol8 = ncol / 8;
  //
  // then allocate the bit map with the land sea dataflags
  // -----------------------------------------------------
  //
  buffer = calloc(nrow,sizeof(unsigned char*));
  for (i=0; i<nrow; i++) {
    buffer[i] = calloc(ncol8,sizeof(unsigned char));
  }
  //
  // then read the mask file
  // -----------------------
  //
  for (i=0; i<nrow; i++) {
    n = fread(buffer[i],ncol8,1,fp);
    if (n != 1) {
      printf("Your mask file is corrupt!\n");
      exit(-1);
    }
  }
  // 
  fclose(fp);
  //printf("Closing file %s\n",inputfile);
  //
  maskset = calloc(8,sizeof(unsigned char));
  maskset[0] = 0x80;
  maskset[1] = 0x40;
  maskset[2] = 0x20;
  maskset[3] = 0x10;
  maskset[4] = 0x08;
  maskset[5] = 0x04;
  maskset[6] = 0x02;
  maskset[7] = 0x01;
  //
  bpd = ncol / 360.0;
}
//
//------------------------------------------------------------------------------
//
void bitmaskend() 
{
  int4 i;
  free(maskset);
  for (i=0; i<nrow; i++) {
    free(buffer[nrow-i-1]);
  }
  free(buffer);
  ncol = 0;
}
//
//------------------------------------------------------------------------------
//
int4 bitmaskeval(real8 *lat, real8 *lon) 
{
  real8 latitude, longitude, colat;
  int4  ix, iy, iz;
  //
  if (ncol == 0) {
    printf("Use bitmaskinit before bitmaskeval\n");
    exit(-1);
  }
  //
  latitude = *lat;
  colat = 90.0-latitude;
  iy = floor(colat * bpd); 
  if ((iy < 0) || (iy >= nrow)) { return -1; }
  //
  longitude = *lon;
  ix = floor(longitude * bpd);
  iz = ix % ncol;
  if (iz < 0) { iz = iz+ncol; }
  //
  if ((buffer[iy][iz/8] & maskset[iz % 8])) {
    return 1;
  } else {
    return 0;
  }
}
//
//------------------------------------------------------------------------------
//
