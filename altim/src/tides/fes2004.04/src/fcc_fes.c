/* ######################################################################
 *  C program for F77/C interface for the FES prediction software.
 * 
 *  File      : fcc_fec.c
 *  Developer : CLS
 *  Version   : 1.1
 *  Date      : 6 oct 2004
 *  
 * ######################################################################
 */       

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "fes.h"

/*
// ///////////////////////////////////////////////////////////////////////////
// 			P R I V A T E  D A T A
// ///////////////////////////////////////////////////////////////////////////
*/

/* Globals variables */
static fesData*	shortTide;
static fesData*	radialTide;

/*
// ///////////////////////////////////////////////////////////////////////////
// 			P R I V A T E  F U N C T I O N S
// ///////////////////////////////////////////////////////////////////////////
*/

static char* rtrim (char* s)
{
  char *p1 = s;

  while (*p1 && !isspace(*p1))
    p1++;

  *p1 = '\0';

  return s;
}

/*
// ///////////////////////////////////////////////////////////////////////////
// Driver to allocate memory used by FES grids
//
// Parameters:
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.
 */
static int initFes(const char* const dir, const tideModel model)
{
  int rc;
  
  rc = fesNew(&shortTide, FES_SHORT_TIDE, model, dir);
  if ( rc != FES_SUCCESS )
    return rc;

  rc = fesNew(&radialTide, FES_RADIAL_TIDE, model, dir);
  return rc;
}


/*
// ///////////////////////////////////////////////////////////////////////////
// Driver to deallocate memory used by FES grids
*/
static void freeFes(void)
{
  fesDelete(shortTide);
  fesDelete(radialTide);
}



/*
// ///////////////////////////////////////////////////////////////////////////
// Driver for tidal computation.
//
// Parameters:
//   handle:	        handle of the internal data->
//   lat		Latitude in degrees (positive north) for the position
//			at which tide is computed.
//   lon		Longitude in degrees (positive north) for the position
//			at which tide is computed.
//   time		Julian day.
//   tide		Computed tide or load tide, in centimeters.
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.
*/
static int fesTide(const double lat,
	    const double lon,
	    const double time,
	    double* tide,
            double* load,
            double* lp)
{
  double loadlp;
  int rc;
  
  if(!shortTide || !radialTide)
  {
    printf("Workspace not initialized.\n");
    return FES_ACCESS_ERROR;
  }

  /* Compute tide */
  rc = fesCore(shortTide, lat, lon, time, tide, lp);
  if ( rc == FES_SUCCESS )
  {
    /* Compute load tide */
    rc = fesCore(radialTide, lat, lon, time, load, &loadlp);
  }
  return rc;
}






/*
// ///////////////////////////////////////////////////////////////////////////
// 			I N T E R F A C E    F O R T R A N  - >  C
// ///////////////////////////////////////////////////////////////////////////
*/


void fcc_initfes_(char* dir, int* model, int* status, int dirLen);
void fcc_initfes_(char* dir, int* model, int* status, int dirLen)
{
  tideModel m	= (tideModel) *model;
  char*	    s	= strdup(dir);
  
  if( s == NULL )
    *status = FES_NO_MEMORY;
  else
  {
    
    *status = initFes(rtrim(s), m);
    free(s);
  }
}

void fcc_festide_(double* lat, 
		  double* lon, 
		  double* time,
	    	  double* tide,
	    	  double* load,
	    	  double* lp,
		  int*    status);
void fcc_festide_(double* lat, 
		  double* lon, 
		  double* time,
	    	  double* tide,
	    	  double* load,
	    	  double* lp,
		  int*    status)
{
/*
 *   printf ("1: %f %f %f\n", *lat, *lon, *time);
 */
  *status = fesTide(*lat, *lon, *time, tide, load, lp);
/*
 *   printf ("2: %f %f %f\n", *tide, *load, *lp);
 */
}

void fcc_freefes_(void);
void fcc_freefes_(void)
{
  freeFes();
}

/* eof festide.c */
