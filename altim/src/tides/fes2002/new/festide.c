/* ######################################################################
 *  Drivers for the FES prediction software.
 * 
 *  Version 6.1  : June, 4 2003
 *                 New version in C (thanks Frédéric Briol !)
 *                 Interface in F77
 *  Major changes: revision by Fabien LEFEVRE and Frédéric BRIOL
 * 
 *  Version 5.2  : May, 21 2003
 *                 Bug with M11 and M12 waves (thanks Shailen.D Desai !)
 *  Major changes: revision by Fabien LEFEVRE
 * 
 *  Version 5.1  : May, 6 2003
 *                 Bug with loading effects (wrong shift in longitude)
 *                 Bug with interpolation
 *  Major changes: revision by Fabien LEFEVRE and Frederic BRIOL
 * 
 *  Version 4.1  : February, 13 2002
 *                 FES2002 version
 *                 Take into account the P1 grid solution
 *  Major changes: revision by Fabien LEFEVRE
 * 
 *  Version 3.2  : May, 9 2001
 *                 update real*4, real*8, integer*4 to be compatible
 *                 with most of Unix and Linux platform 
 *  Minor changes: revision by Fabien LEFEVRE
 * 
 *  Version 3.1  : November, 1 2000
 *  Major changes: full revision by Fabien LEFEVRE
 * 
 * 
 *  This software have been tested on SUN Solaris 2.7
 *  and Linux platform
 *  It is provided without any guarantees.
 * 
 *  For bug reports, please contact :
 *  ---------------------------------
 *  Fabien LEFEVRE 
 * 
 *  CLS
 *  http://www.cls.fr
 *  Direction Océanographie Spatiale
 *  8-10, rue Hermès - Parc Technologique du Canal
 *  31526 Ramonville Saint-Agne cedex - France
 *  Tel: +33 (0)5 61 39 37 45 Fax: +33 (0)5 61 39 37 82
 *  e-mail: Fabien.Lefevre@cls.fr
 * 
 *  NOTE: This software is based on the former versions
 *        developed by Jean-Marc MOLINES and Florent LYARD
 * #####################################################################*/

#include <assert.h>
#include <stdio.h>

#include "festide.h"

/*
// ///////////////////////////////////////////////////////////////////////////
// 			P R I V A T E  D A T A
// ///////////////////////////////////////////////////////////////////////////
*/

/* Globals variables */
static void*	hTide;
static void*	hLoad;
#define		TIDE	0
#define		LOAD	1
#define		VERBOSE	1

/*
// ///////////////////////////////////////////////////////////////////////////
// 			P U B L I C S  F U N C T I O N S
// ///////////////////////////////////////////////////////////////////////////
*/

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
//   1 is returned to indicate an error otherwise 0.
*/
int initFes(const char* const dir, const tideModel model)
{
  if(newHandle(&hTide, fes2002, TIDE, VERBOSE))
    goto onError;

  if(newHandle(&hLoad, fes2002, LOAD, VERBOSE))
    goto onError;
  
  if(loadGrids(hTide, dir))
    goto onError;

  if(loadGrids(hLoad, dir))
    goto onError;

  return 0;
  
onError:
  return 1;
}


/*
// ///////////////////////////////////////////////////////////////////////////
// Driver to deallocate memory used by FES grids
//
// Parameters:
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
void freeFes(void)
{
  freeHandle(hTide);
  freeHandle(hLoad);
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
//   1 is returned to indicate an error otherwise 0.
*/
int fesTide(const double lat,
	    const double lon,
	    const double time,
	    double* tide,
            double* load,
            double* lp)
{
  if(!hTide || !hLoad)
  {
    fprintf(stderr, "Workspace not initialized.\n");
    goto onError;
  }

  if(fesCore(hTide, lat, lon, time, tide))
    goto onError;

  if(fesCore(hLoad, lat, lon, time, load))
    goto onError;

  lpeqmt(time, lat, lp);

  return 0;
onError:
  return 1;
}






/*
// ///////////////////////////////////////////////////////////////////////////
// 			I N T E R F A C E    F O R T R A N  - >  C
// ///////////////////////////////////////////////////////////////////////////
*/

#ifdef NOUNDERSCOREAFTER
void fccinitfes(char* dir, int* model, int* status, int dirLen);
void fccinitfes(char* dir, int* model, int* status, int dirLen)
#else
void fccinitfes_(char* dir, int* model, int* status, int dirLen);
void fccinitfes_(char* dir, int* model, int* status, int dirLen)
#endif
{
  tideModel m	= (tideModel) *model;
  *status = initFes(rtrim(dir), m);
  return;
}

#ifdef NOUNDERSCOREAFTER
void fccfestide(double* lat, 
		  double* lon, 
		  double* time,
	    	  double* tide,
	    	  double* load,
	    	  double* lp,
		  int*    status);
void fccfestide(double* lat, 
		  double* lon, 
		  double* time,
	    	  double* tide,
	    	  double* load,
	    	  double* lp,
		  int*    status)
#else
void fccfestide_(double* lat, 
		  double* lon, 
		  double* time,
	    	  double* tide,
	    	  double* load,
	    	  double* lp,
		  int*    status);
void fccfestide_(double* lat, 
		  double* lon, 
		  double* time,
	    	  double* tide,
	    	  double* load,
	    	  double* lp,
		  int*    status)
#endif
{
  *status = fesTide(*lat, *lon, *time, tide, load, lp);
  return;
}

#ifdef NOUNDERSCOREAFTER
void fccfreefes(void);
void fccfreefes(void)
#else
void fccfreefes_(void);
void fccfreefes_(void)
#endif
{
  freeFes();
  return;
}

/* eof festide.c */
