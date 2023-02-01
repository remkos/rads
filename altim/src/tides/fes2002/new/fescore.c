/* ######################################################################
 *  Main routines for the FES prediction software.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "fescore.h"

/*
// ///////////////////////////////////////////////////////////////////////////
// 			P R I V A T E  D A T A
// ///////////////////////////////////////////////////////////////////////////
*/

#define DV		18446744073709551616.0
#define	MAX_PATH  	1024
#define IN_GAP(a, b, c)	(((b) >= (a)) && ((b) <= (c)))
#define DEPTH_FES2002	18
#define DEPTH_FES99	16
#define EPSILON		1.0E-9
#define NB_WAVE		27
#define DEPTH		DEPTH_FES2002
#define SQR(a)		((a)*(a))
#define MIN(a, b)	((a) < (b)? (a): (b))
#define EQUALS(a, b)	(fabs(a - b) < EPSILON)

/*
// ///////////////////////////////////////////////////////////////////////////
// Definition of a grids
*/
typedef struct
{
  double 	latMin;
  double	latMax;
  double 	lonMin;
  double	lonMax;
  double 	latStep;
  double 	lonStep;
  int	 	latSamples;
  int	 	lonSamples;
  int		depth;
  float***	data;
} gridDsc;

/*
// ///////////////////////////////////////////////////////////////////////////
// Internal data
*/
typedef struct
{
  short		inGrid;
  short		isData;
  short		computeP1;
  short		shortTide;
  short		verbose;
  short		num[NB_WAVE];
  double	aamu2;
  double	aanu2;
  double	aal2;
  double	aat2;
  double	aalda2;
  double	aap1;
  double	bbmu2;
  double	bbnu2;
  double	bbl2;
  double	bbt2;
  double	bblda2;
  double	bbp1;
  double	ccmu2;
  double	ccnu2;
  double	ccl2;
  double	cct2;
  double	cclda2;
  double	ccp1;
  double	tNodal;
  double	southLat;
  double	northLat;
  double	westLon;
  double	eastLon;
  double	sw[DEPTH];
  double	se[DEPTH];
  double	nw[DEPTH];
  double	ne[DEPTH];
  double	freq[NB_WAVE];
  double	f[NB_WAVE];
  double	v0u[NB_WAVE];
  double	re[NB_WAVE];
  double	im[NB_WAVE];
  gridDsc	grid;
} fesData;

static const int quantiemes[2][12] =
{
  { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 },
  { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 }
};

static const char* waveName[] = { "Q1", "O1", "K1", "2N2", "N2", "M2", "K2", "S2", "P1", NULL };

/* Globals constants */
static double	pi	= 3.14159265358979323846;
static double	rad	= 3.14159265358979323846 / 180.0;
static double	deg	= 180.0 / 3.14159265358979323846;

#ifdef WIN32
/* NTFS file system */
static const char* tideFmt	= "%s\\%s_%s.asc";
static const char* loadFmt	= "%s\\%s_dr%s.asc";
#else
/* UNIX file system */
static const char* tideFmt	= "%s/%s_%s.asc";
static const char* loadFmt	= "%s/%s_dr%s.asc";
#endif





/*
// ///////////////////////////////////////////////////////////////////////////
// 			P R I V A T E S  F U N C T I O N S
// ///////////////////////////////////////////////////////////////////////////
*/






/*
// ///////////////////////////////////////////////////////////////////////////
// Standardization of longitude
//
// Parameters:
//	base:	base standardization.
//	lon:	longitude
//
// Return value:
//   Standardized longitude.
*/
static double normalizeLongitude(const double base, const double lon)
{
  register double result = lon;

  while(result >= ((base + 360.0) - EPSILON))
    result -= 360.0;

  while(result < base - EPSILON)
    result += 360.0;

  if(fabs(result - base) <= EPSILON)
    result = base;

  return result;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Whole value nearest to a real.
//
// Parameters:
//	x:	real
//
// Return value:
//   Whole value
*/
static int round(const double x)
{
  double  i;
  double  f = modf(x, &i);

  return (int) (fabs(f) > 0.5? i + (f < 0.0? -1.0: 1.0): i);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Fractional part of a real.
//
// Parameters:
//	x:	real
//
// Return value:
//   Fractional part
*/
static double frac(const double x)
{
  double i;

  return modf(x, &i);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Linear Weighting (M5)
//
// Parameters:
//	x:	abscissa where interpolation will occur
//	x1: 	abscissa corresponding to the first value
//	x2:	abscissa corresponding to the second value
//	w1	linear weight w1
//	w2:	linear weight w2
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
static int linearWeighting(const double x,
			   const double x_1,
			   const double x_2,
			   double* w_1,
			   double* w_2)
{
  if((x < x_1 - EPSILON) || (x > x_2 + EPSILON))
  {
    fprintf(stderr, "Coordinate (%f) out of bounds [%f,%f]\n", x, x_1, x_2);
    return 1;
  }

  if(EQUALS(x_1, x_2) || EQUALS(x, x_1))
  {
    *w_1 = 1.0;
    *w_2 = 0.0;
  }
  else if (EQUALS(x, x_2))
  {
    *w_1 = 0.0;
    *w_2 = 1.0;
  }
  else
  {
    *w_1 = (x_2 - x) / (x_2 - x_1);
    *w_2 = (x - x_1) / (x_2 - x_1);
  }

  return 0;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// bilinear interpolation.
//
// Parameters:
//	x_1:		X-coordinate X1
//	x_2: 		X-coordinate X2
//	y_1:		Y-coordinate Y1
//	y_2		Y-coordinate Y2
//	value_11:	Value of the point (X1, Y1)
//	value_21:	Value of the point (X2, Y1)
//	value_12:	Value of the point (X1, Y2)
//	value_22:	Value of the point (X2, Y2)
//	x:		X-coordinate of the point where the interpolation is
//			carried out
//	x:		Y-coordinate of the point where the interpolation is
//			carried out
//	result:		Result of the interpolation.
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
static int bilinearInterp(const double x_1,
			  const double x_2,
			  const double y_1,
			  const double y_2,
			  const double value_11,
			  const double value_21,
			  const double value_12,
			  const double value_22,
			  const double x,
			  const double y,
			  double* result)
{
  double	w_X1;
  double	w_X2;
  double	w_Y1;
  double	w_Y2;
  double	w_11, w_21, w_12, w_22;
  double	w				= 0.0;

  if (linearWeighting(x, x_1, x_2, &w_X1, &w_X2))
    goto onError;

  if (linearWeighting(y, y_1, y_2, &w_Y1, &w_Y2))
    goto onError;

  if (EQUALS(value_11, DV))
  {
    w_11 = 0.0;
  }
  else
  {
    w_11 = w_X1 * w_Y1;
    w   += w_11;
  }

  if (EQUALS(value_21, DV))
  {
    w_21 = 0.0;
  }
  else
  {
    w_21 = w_X2 * w_Y1;
    w   += w_21;
  }

  if (EQUALS(value_12, DV))
  {
    w_12 = 0.0;
  }
  else
  {
    w_12 = w_X1 * w_Y2;
    w   += w_12;
  }

  if (EQUALS(value_22, DV))
  {
    w_22 = 0.0;
  }
  else
  {
    w_22 = w_X2 * w_Y2;
    w   += w_22;
  }

  *result = (w_11 * value_11 + w_21 * value_21 + w_12 * value_12 + w_22 * value_22) / w;

  return 0;

onError:
  return 1;
}










/*
// ///////////////////////////////////////////////////////////////////////////
// Whole part of a real
//
// Parameters:
//	x:	real
//
// Return value:
//   Whole part.
*/
static double nInt(const double x)
{
  double integral;

  modf(x, &integral);

  return integral;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Conversion of Seconds 1950 into calendar date..
//
// Parameters:
//	date:		Number of seconds passed since 1950.
//	day:		Day in the month.
//	month:		Month in the year.
//	year:		Year (4 digits)
//	hours:		hours in the day
//	minuts:		Minuts
//	seconds:	Seconds
//	mcsec:		micro seconds
*/
static void decodeDate(const double date,
		       int* day,
		       int* month,
		       int* year,
		       int* hours,
		       int* minuts,
		       int* seconds,
		       int* mcsec)
{
  double	tmp;
  int		nDays;
  int		nYears;
  int		nYearsBis;
  int		iYearsBis;

  assert(date >= 0);

  *mcsec = round(frac(date) * 1000000.0);

  tmp = nInt(date);

  nDays = (int) (tmp / 86400.0);
  tmp  -= nDays * 86400.0;

  *hours = (int) nInt(tmp / 3600.0);
  tmp -= *hours * 3600.0;

  *minuts = (int) nInt(tmp / 60.0);
  tmp -= *minuts * 60.0;

  *seconds = (int) nInt(tmp);

  if (*mcsec >= 1000000)
  {
    *mcsec -= 1000000;
    if (++(*seconds) > 60)
    {
      *seconds -= 60;
      if (++(*minuts) >= 60)
      {
	*minuts	-= 60;
	if (++(*hours) >= 24)
	{
	  *hours -= 24;
	  nDays++;
	}
      }
    }
  }

  tmp = nDays;
  nYears = (int) nInt((tmp + 0.5) / 365.25);

  *year	= nYears + 1950;

  nYearsBis =  (*year - 1 - 1900) / 4 - 12;

  iYearsBis = MIN((*year % 4) + 1, 2) - 1;

  nDays	= nDays - (nYears * 365) - nYearsBis + 1;

  *month = 1;

  while ((*month < 12) && (quantiemes[iYearsBis][*month] < nDays))
    (*month)++;

  *day	= nDays -  quantiemes[iYearsBis][*month - 1];
}




/*
// ///////////////////////////////////////////////////////////////////////////
// Returns the index in the grid corresponding to a value
//
// Parameters:
//	min:	Min value in the grid
//	value:  Asked value
//	step:	Step
//
// Return value:
//   Computed index.
*/
static int getIndex(const double min, const double value, const double step)
{
  return (int) nInt((value - min) / step);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Returns a value corresponding to an index
//
// Parameters:
//	idx:	Asked index
//	min:	Min value in the grid
//	step:	Step
//
// Return value:
//   Computed value.
*/
static double getValue(const int idx, const double min, const double step)
{
  return min + (step * idx);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Allocates a grid in memory with elements initialized to 0.
//
// Parameters:
//   grid:   grid to allocate
//
// Return value:
//   1 is returned to indicate an error, otherwise 0.
*/
static int newGrid(gridDsc* grid)
{
  int i, j;

  if((grid->data = (float***) calloc(grid->lonSamples, sizeof(float))) == NULL)
    goto onError;

  for(i = 0; i < grid->lonSamples; ++i)
  {
    if((grid->data[i] = (float**) calloc(grid->latSamples, sizeof(float))) == NULL)
      goto onError;

    for(j = 0; j < grid->latSamples; ++j)
    {
      if((grid->data[i][j] = (float*) calloc(grid->depth, sizeof(float))) == NULL)
	goto onError;
    }
  }

  return 0;

onError:
  fprintf(stderr,
    "Not enought memory to allocate %d bytes.\n",
    grid->lonSamples * grid->latSamples * grid->depth * sizeof(float));
  return 1;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Read formatted data from a line.
//
// Parameters:
//   line:	line to read
//   undef:	Value of a cell nondefined in the grid.
//   x:		longitude to read
//   y:		latitude to read
//   depth:	Current grid depth
//   grid:	Grid containing the data read.
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
static int readValues(char* line,
		      const float undef,
		      int* x,
		      const int y,
		      const int depth,
		      gridDsc* grid)
{
  char	*token = strtok(line, " \t");

  while(token)
  {
    if(sscanf(token, "%f", &grid->data[*x][y][depth]) != 1)
      return 1;

    if(grid->data[*x][y][depth] == undef)
      grid->data[*x][y][depth] = DV;

    (*x)++;

    token  = strtok(NULL, " \t");
  }

  return 0;
}






/*
// ///////////////////////////////////////////////////////////////////////////
// Read formatted grid from filename.
//
// Parameters:
//   filename	File containing the grid.
//   depth	current depth
//   model:	1 reading fes99 ortherwise fes2002
//   grid:	Grid containing the data read.
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
static int readGrid(const char* const filename,
		    const int depth,
		    const short model,
		    gridDsc* grid)
{
  FILE*		file;
  char		buffer[MAX_PATH];
  int		result 			= 0;
  int		nX;
  int		nY;
  int		x;
  int		y			= 0;
  float		undefX;
  float		undefY;
  double	unused;
  double	latMin;
  double	lonMin;
  double	latStep;
  double	lonStep;

  if((file = fopen(filename, "rt")) == NULL)
  {
    fprintf(stderr, "Can't read %s.\n", filename);
    goto onError;
  }

  /* Reading XMIN & XMAX (X <=> Longitude) */
  if(fgets(buffer, MAX_PATH, file) == NULL)
    goto onIOError;

  if(sscanf(buffer, "%lf %lf",&lonMin, &unused) !=  2)
    goto onIOError;

  /* Reading YMIN & YMAY (Y <=> Latitude) */
  if(fgets(buffer, MAX_PATH, file) == NULL)
    goto onIOError;

  if(sscanf(buffer, "%lf %lf",&latMin, &unused) !=  2)
    goto onIOError;

  /* Reading DX & DY */
  if(fgets(buffer, MAX_PATH, file) == NULL)
    goto onIOError;

  if(sscanf(buffer, "%lf %lf",&lonStep, &latStep) !=  2)
    goto onIOError;

  /* Reading NX & DY */
  if(fgets(buffer, MAX_PATH, file) == NULL)
    goto onIOError;

  if(sscanf(buffer, "%d %d",&nX, &nY) !=  2)
    goto onIOError;

  /* The grids should not overlap. */
  nX--;
  nY--;

  /* Reading MASKX, MASKY */
  if(fgets(buffer, MAX_PATH, file) == NULL)
    goto onIOError;

  if(sscanf(buffer, "%f %f",&undefX, &undefY) != 2)
    goto onIOError;

  assert(undefX == undefY);

  /* First call */
  if(depth == 0)
  {
    grid->latMin	= latMin;
    grid->lonMin	= lonMin;
    grid->latStep	= latStep;
    grid->lonStep	= lonStep;
    grid->latSamples	= nY;
    grid->lonSamples	= nX;
    grid->depth		= model? DEPTH_FES99: DEPTH_FES2002;

    if(newGrid(grid))
      goto onError;
  }
  /* Check grid definition */
  else if(grid->latMin	!= latMin	|| grid->lonMin		!= lonMin ||
    grid->latStep	!= latStep	|| grid->lonStep	!= lonStep ||
    grid->latSamples	!= nY		|| grid->lonSamples	!= nX)
  {
    fprintf(stderr, "The Definition of the grids is not constant in the input files.\n");
    fprintf(stderr, "Awaited definition:\n");
    fprintf(stderr, "\tlatMin: %f\tlonMin: %f\n", latMin, lonMin);
    fprintf(stderr, "\tlatStep: %f\tlonStep: %f\n", latStep, lonStep);
    fprintf(stderr, "\tlatSamples: %d\t\tlonSamples: %d\n", nY, nX);
    fprintf(stderr, "Definition read:\n");
    fprintf(stderr, "\tlatMin: %f\tlonMin: %f\n", grid->latMin, grid->lonMin);
    fprintf(stderr, "\tlatStep: %f\tlonStep: %f\n", grid->latStep, grid->lonStep);
    fprintf(stderr, "\tlatSamples: %d\t\tlonSamples: %d\n", grid->latSamples, grid->lonSamples);
    goto onError;
  }

  /* For all the latitudes. */
  while(y < nY)
  {
    x = 0;

    /* For all the longitudes. */
    while(x < nX)
    {
      /* Reading the 30 values of the amplitude. */
      if(fgets(buffer, MAX_PATH, file) == NULL)
	goto onIOError;

      if(readValues(buffer, undefY, &x, y, depth, grid))
	goto onIOError;

      x -= 30;

      /* Reading the 30 values of the phase. */
      if(fgets(buffer, MAX_PATH, file) == NULL)
	goto onIOError;

      if(readValues(buffer, undefY, &x, y, depth + 1, grid))
	goto onIOError;
    }
    /* Reading of the redundant value for the amplitude. */
    if(fgets(buffer, MAX_PATH, file) == NULL)
      goto onIOError;

    /* Reading of the redundant value for the phase. */
    if(fgets(buffer, MAX_PATH, file) == NULL)
      goto onIOError;

    y++;
  }

  goto onTerminate;

onIOError:
  fprintf(stderr, "Error reading file %s.\n", filename);

onError:
  result = 1;

onTerminate:
  if(file)
    fclose(file);
  return result;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// This program initialize some astronomic data useful for
// nodal corrections.
//   itj	Desired UTC time, in (decimal) Modified.
//   tt		Mean solar angle relative to Greenwich
//   hp		Mean solar longitude
//   s		Mean lunar longitude
//   p1		Longitude of solar perigee
//   p		Longitude of lunar perigee</param>
//   iang
//   xi
//   nu
//   x1ra
//   r
//   nuprim
//   nusec
*/
static void astronomics(const double  itj,
			double*	tt,
			double*	hp,
			double*	s,
			double*	p1,
			double*	p,
			double*	iang,
			double*	xi,
			double*	nu,
			double*	x1ra,
			double*	r,
			double*	nuprim,
			double*	nusec)
{
  static const double	ct0 	     = 180.0;
  static const double	ct1 	     = 360.0 * 3.6525E+04;
  static const double	cn0 	     = 259.1560563;
  static const double	cn1 	     = -1934.1423972;
  static const double	cs0 	     = 277.0256206;
  static const double	cs1 	     = 481267.892;
  static const double	ch0 	     = 280.1895015;
  static const double	ch1 	     = 36000.76892;
  static const double	cps0	     = 281.2208568;
  static const double	cps1	     = 1.719175;
  static const double	cp0 	     = 334.3837214;
  static const double	cp1 	     = 4069.0322056;
  double		tgn2;
  double		at1;
  double		at2;
  double		u;
  double		tgi2;
  double		n;
  double		pp;

  /* tt mean solar angle relative to Greenwich */
  *tt	= fmod(ct0 + ct1 * itj , 360.0);

  /* hp longitude of ascending lunar node */
  n	= fmod(cn0 + cn1 * itj, 360.0) * rad;

  /* hp mean solar longitude */
  *hp	= fmod(ch0 + ch1 * itj, 360.0) * rad;

  /* s mean lunar longitude */
  *s	= fmod(cs0 + cs1 * itj, 360.0) * rad;

  /* p1 longitude of solar perigee */
  *p1	= fmod(cps0 + cps1 * itj, 360.0) * rad;

  /* p longitude of lunar perigee */
  *p	= fmod(cp0 + cp1 * itj, 360.0) * rad;

  u	= 9.13694997e-01 - 3.5692561e-02 * cos(n);

  *iang	= acos(u);

  tgn2	= tan(n / 2.0);

  at1	= atan(1.01883 * tgn2);
  at2	= atan(6.4412e-01 * tgn2);

  *xi = -at1 - at2 + n;
  if(n > pi)
  {
    *xi -= 2.0 * pi;
  }

  *nu = at1 - at2;

  /* for constituents l2,k1,k2 */
  tgi2	= tan(*iang / 2.0);

  pp	= *p - *xi;

  *x1ra	= sqrt(1.0 - 12.0 * SQR(tgi2) * cos(2.0 * pp) + 36.0 * SQR(SQR(tgi2)));

  *r	= atan(sin(2.0 * pp) / (1.0/(6.0 * SQR(tgi2)) - cos(2.0 * pp)));

  *nuprim =
    atan(sin(2.0 * (*iang)) * sin(*nu) /
    (sin(2.0 * (*iang)) * cos(*nu) + 3.347E-01));

  *nusec = 0.5 *
    atan(((SQR(sin(*iang))) * sin(2.0 * (*nu))) /
    (SQR(sin(*iang)) * cos(2.0 * (*nu))+ 7.27E-02));
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Compute nodal corrections from SCHUREMAN (1958).
// indexes used in this routine are internal to the code
// and corresponds to the "original" ondes.dat file.
//
// Parameters:
//   data:	Internal data->
//   iang
//   nu
//   x1ra
*/
static void nodalA(fesData* data,
		   const double	iang,
		   const double nu,
		   const double x1ra)
{
  int	i;

  for(i = 0; i < NB_WAVE; i++)
  {
    switch(data->num[i])
    {
    case 1:
    case 27:
    case 65:
    case 66:
    case 67:
    case 69:
      data->f[i] = sin(iang) * SQR(cos(iang / 2.0))/0.38;
      break;
    case 2:
    case 12:
    case 13:
    case 60:
    case 71:
    case 72:
      data->f[i] = 1.0;
      break;
    case 3:
      data->f[i] =
	sqrt(0.8965 * SQR(sin(2.0 * iang)) +
	0.6001 * sin(2.0 * iang) * cos(nu) + 0.1006);
      break;
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 61:
    case 62:
      data->f[i] = SQR(SQR(cos(iang / 2.0))) / 0.9154;
      break;
    case 11:
      data->f[i] = SQR(SQR(cos(iang / 2.0))) / 0.9154 * x1ra;
      break;
    case 14:
      data->f[i] = sqrt(19.0444 * SQR(SQR(sin(iang))) +
	2.7702 * SQR(sin(iang)) * cos(2.0 * nu) + 0.0981);
      break;
    case 63:
    case 64:
      data->f[i] = sin(iang) * sin(iang) / 0.1565;
      break;
    case 68:
    case 70:
    case 73:
    case 74:
      data->f[i] = sin(2.0 * iang) / 0.7214;
      break;
    case 75:
      data->f[i] = sin(iang) * SQR(sin(iang / 2.0)) / 0.01640;
      break;
    default:
      assert(0);
    }
  }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Compute nodal corrections from SCHUREMAN (1958).
// indexes used in this routine are internal to the code
// and corresponds to the "original" ondes.dat file.
//
// Parameters:
//  tt
//  hp
//  s
//  p1
//  p
//  xi
//  nu
//  r
//  nupri
//  nusec
*/
static void nodalG(fesData* data,
		   const double	tt,
		   const double	hp,
		   const double	s,
		   const double	p1,
		   const double	p,
		   const double	xi,
		   const double	nu,
		   const double	r,
		   const double	nuprim,
		   const double	nusec)
{
  int		i;
  double	u;
  double	v0;

  for(i = 0; i < NB_WAVE; i++)
  {
    switch(data->num[i])
    {
      /* O1 */
    case 1:
      /*v0 = tt - 2.0 * s + hp + 90.0;
      u = 2 * xi - nu;*/
      data->v0u[i] = tt - 2.0 * s + hp+ 90.0 +
	2.0 * xi - nu;
      break;
      /* P1 */
    case 2:
      /*v0 = tt - hp + 90.0;
      u = 0.0;*/
      data->v0u[i]=tt - hp + 90.0;
      break;
      /* K1 */
    case 3:
      /*v0 = tt + hp - 90.0;
      u = -nuprim;*/
      data->v0u[i] = tt + hp - 90.0 - nuprim;
      break;
      /* 2N2 */
    case 5:
      /*v0 = 2.0 * tt - 4.0 * s + 2.0 * hp + 2.0 * p;
      u = 2.0 * xi - 2.0 * nu;*/
      data->v0u[i] =2.0 * tt - 4.0 * s + 2.0 * hp +
	2.0 * p + 2.0 * xi - 2.0 * nu;
      break;
      /* MU2 */
    case 6:
      /*v0 = 2.0 * tt - 4.0 * s + 4.0 * hp;
      u = 2.0 * xi - 2.0 * nu;*/
      data->v0u[i] = 2.0 * tt - 4.0 * s + 4.0 * hp +
	2.0 * xi - 2.0 * nu;
      break;
      /* N2 */
    case 7:
      /*v0 = 2.0 * tt - 3.0 * s + 2.0 * hp + p;
      u = 2.0 * xi - 2.0 * nu;*/
      data->v0u[i] = 2.0 * tt - 3.0 * s + 2.0 * hp +
	p + 2.0 * xi - 2.0 * nu;
      break;
      /* Nu2 */
    case 8:
      /*v0 = 2.0 * tt - 3.0 * s + 4.0 * hp - p;
      u = 2.0 * xi - 2.0 * nu;*/
      data->v0u[i] = 2.0 * tt - 3.0 * s + 4.0 * hp -
	p + 2.0 * xi - 2.0 * nu;
      break;
      /* M2 */
    case 9:
      /*v0 = 2.0 * tt - 2.0 * s + 2.0 * hp;
      u = 2.0 * xi - 2.0 * nu;*/
      data->v0u[i] = 2.0 * tt + - 2.0 * s + 2.0 * hp +
	2.0 * xi - 2.0 * nu;
      break;
      /* L2 */
    case 11:
      /*v0 = 2.0 * tt - s + 2.0 * hp - p + 180.0;
      u = 2.0 * xi - 2.0 * nu - r;*/
      data->v0u[i] = 2.0 * tt - s + 2.0 * hp -
	p + 180.0 + 2.0 * xi - 2.0 * nu - r;
      break;
      /* T2 */
    case 12:
      /*v0 = 2.0 * tt - hp + p1;
      u = 0.0;*/
      data->v0u[i] = 2.0 * tt - hp + p1;
      break;
      /* S2 */
    case 13:
      /*v0 = 2.0 * tt;
      u = 0.0;*/
      data->v0u[i] = 2.0 * tt;
      break;
      /* K2 */
    case 14:
      /*v0 = 2.0 * tt + 2.0 * hp;
      u = -2.0 * nusec;*/
      data->v0u[i] = 2.0 * tt + 2.0 * hp - 2.0 * nusec;
      break;
      /* Q1 */
    case 27:
      /*v0 = tt - 3.0 * s + hp + p + 90.0;
      u = 2.0 * xi - nu;*/
      data->v0u[i] = tt - 3.0 * s + hp + p + 90.0 + 2.0 * xi - nu;
      break;
      /* Eps2 */
    case 60:
      /*v0 = 2.0 * tt - 5 * s + 4 * hp + p;
      u = 0.0;*/
      data->v0u[i] = 2.0 * tt - 5.0 * s + 4.0 * hp + p;
      break;
      /* Lambda2 */
    case 61:
      v0 = 2.0 * tt - s + p + 180.0;
      u = 2.0 * xi - 2.0 * nu;
      data->v0u[i] = v0 + u;
      break;
      /* l21 */
    case 62:
      v0 = 2.0 * tt - s + 2.0 * hp - p + 180.0;
      u = -2.0 * nu;
      data->v0u[i] = v0 + u;
      break;
      /* l22 */
    case 63:
      v0 = 2.0 * tt - s + 2.0 * hp + p;
      u = 2.0 * xi - 2.0 * nu;
      data->v0u[i] = v0 + u;
      break;
      /* Eta2 */
    case 64:
      v0 = 2.0 * tt + s + 2.0 * hp - p;
      u = -2.0 * nu;
      data->v0u[i] = v0 + u;
      break;
      /* 2Q1 */
    case 65:
      v0 = tt - 4.0 * s + hp + 2.0 * p + 90.0;
      u = 2.0 * xi - nu;
      data->v0u[i] = v0 + u;
      break;
      /* Sigma1 */
    case 66:
      v0 = tt - 4.0 * s + 3.0 * hp + 90.0;
      u = 2.0 * xi - nu;
      data->v0u[i] = v0 + u;
      break;
      /* Ro1 */
    case 67:
      v0 = tt - 3.0 * s + 3.0 * hp - p + 90.0;
      u = 2.0 * xi - nu;
      data->v0u[i] = v0 + u;
      break;
      /* M11 */
    case 68:
      v0 = tt - s + hp + p - 90.0;
      u = -nu;
      data->v0u[i] = v0 + u;
      break;
      /* M12 */
    case 69:
      v0 = tt - s + hp - p - 90.0;
      u = 2.0 * xi - nu;
      data->v0u[i] = v0 + u;
      break;
      /* Ki1 */
    case 70:
      v0 = tt - s + 3.0 * hp - p - 90.0;
      u = -nu;
      data->v0u[i] = v0 + u;
      break;
      /* Pi1 */
    case 71:
      /*v0 = tt - 2.0 * hp + p1 + 90.0;
      u = 0.0;*/
      data->v0u[i] = tt - 2.0 * hp + p1 + 90.0;
      break;
      /* Phi1 */
    case 72:
      /*v0 = tt + 3.0 * hp - 90.0;
      u = 0.0;*/
      data->v0u[i] = tt + 3.0 * hp - 90.0;
      break;
      /* Teta1 */
    case 73:
      v0 = tt + s - hp + p - 90.0;
      u = -nu;
      data->v0u[i] = v0 + u;
      break;
      /* J1 */
    case 74:
      v0 = tt + s + hp - p - 90.0;
      u = -nu;
      data->v0u[i] = v0 + u;
      break;
      /* OO1 */
    case 75:
      v0 = tt + 2.0 * s + hp - 90.0;
      u = -2.0 * xi - nu;
      data->v0u[i] = v0 + u;
      break;
    default:
      assert(0);
    }
    data->v0u[i] = fmod(data->v0u[i], 360.00);
  }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Nodal correction calls.
//
// Parameters:
//  data	Internal data->
//  tj:		Desired UTC time, in (decimal) Modified.
*/
static void nodalC(fesData* data, const double	tj)
{
  double	tt;
  double	hp;
  double	s;
  double	p1;
  double	p;
  double	iang;
  double	xi;
  double	nu;
  double	x1ra;
  double	r;
  double	nuprim;
  double	nusec;

  astronomics(tj, &tt, &hp, &s, &p1, &p, &iang, &xi, &nu, &x1ra, &r,
    &nuprim, &nusec);

  nodalA(data, iang, nu, x1ra);

  nodalG(data, tt, hp * deg, s * deg, p1 * deg, p * deg, xi * deg, nu * deg,
    r * deg, nuprim * deg, nusec * deg);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Compute nodal corrections.
//
// Parameters:
//  data:	Internal data
//  t0:		Desired UTC time, in (decimal) Modified.
*/
static void initCorrections(fesData* data,
			    const double t0)
{
  int i;

  nodalC(data, t0);

  for(i=0; i < NB_WAVE; i++)
    data->v0u[i] *= rad;

  data->tNodal = t0;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Reading the points closest to a position.
//
// Parameters:
//   grid:	Grid containing the data->
//   lat:	Latitude of the desired position.
//   lon:	Latitude of the desired position.
//   depth:	Depth of the grid.
//   sw:	value for the SW point.
//   se:	value for the SE point.
//   nw:	value for the NW point.
//   ne:	value for the NE point.
//   southLat:  latitude for the southern point.
//   northLat:  latitude for the northern point.
//   westLon:   longitude for the western point.
//   eastLon:   longitude for the eastern point.
//   inGrid:	1 the position is in the grid, otherwise 0.
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
static int getNearestPoints(gridDsc* grid,
			    const double lat,
			    const double lon,
			    const int depth,
			    double* sw,
			    double* se,
			    double* nw,
			    double* ne,
			    double* southLat,
			    double* northLat,
			    double* westLon,
			    double* eastLon,
			    short* inGrid)
{
  double lonNorm = normalizeLongitude(grid->lonMin, lon);

  if(depth != grid->depth)
  {
    fprintf(stderr, "Grid depth is not the expected one.\n"
      "Expected = %d\nCurrent = %d\n", grid->depth, depth);
    return 1;
  }

  /* Check if asked position is in the grid */
  if ((!IN_GAP(grid->latMin, lat, grid->latMax)) ||
      (!IN_GAP(grid->lonMin, lonNorm, grid->lonMax)))
  {
    *inGrid = 0;
  }
  else
  {
    int iLat1 = getIndex(grid->latMin, lat, grid->latStep);
    int iLon1 = getIndex(grid->lonMin, lonNorm, grid->lonStep);
    int iLat2;
    int iLon2;
    int iDepth;

    if(lat >= grid->latMax)
    {
      iLat2 = iLat1;
      iLat1--;
    }
    else
      iLat2 = iLat1 + 1;

    if(lonNorm >= grid->lonMax)
    {
      iLon2 = iLon1;
      iLon1--;
    }
    else
      iLon2 = iLon1 + 1;

    *inGrid   = 1;

    *southLat = getValue(iLat1, grid->latMin, grid->latStep);
    *northLat = getValue(iLat2, grid->latMin, grid->latStep);
    *westLon  = getValue(iLon1, grid->lonMin, grid->lonStep);
    *eastLon  = getValue(iLon2, grid->lonMin, grid->lonStep);

    if(*westLon != lon)
    {
      double gap  = *westLon - *eastLon;
      lonNorm	  = normalizeLongitude(lon, *eastLon);
      *eastLon	  = lonNorm;
      *westLon	  = lonNorm + gap;
    }

    iLon1 %=  grid->lonSamples;
    iLon2 %=  grid->lonSamples;

    for(iDepth = 0; iDepth < depth; iDepth++)
    {
      sw[iDepth] = grid->data[iLon1][iLat1][iDepth];
      se[iDepth] = grid->data[iLon2][iLat1][iDepth];
      nw[iDepth] = grid->data[iLon1][iLat2][iDepth];
      ne[iDepth] = grid->data[iLon2][iLat2][iDepth];
    }
  }
  return 0;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// 			P U B L I C S  F U N C T I O N S
// ///////////////////////////////////////////////////////////////////////////
*/





/*
// ///////////////////////////////////////////////////////////////////////////
// Load grids into memory.
//
// Parameters:
//   handle:	Handle of the internal data.
//   dir:    	directory
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
int loadGrids(void* handle, const char* const dir)
{
  char filename[MAX_PATH];
  int i = 0;
  fesData* data = (fesData*) handle;

  memset(&data->grid, 0, sizeof(gridDsc));

  if(data->verbose)
  {
    fprintf(stdout, "Loading %s grids.\n", data->shortTide?
    	"pure oceanic tide": "radial loading tide");
  }

  while(waveName[i])
  {
    if(!strcmp(waveName[i], "P1") && data->computeP1)
      continue;

    sprintf(filename,
      data->shortTide? tideFmt: loadFmt,
      dir,
      waveName[i],
      data->computeP1? "fes99": "fes2002");

    if(data->verbose)
    {
      fprintf(stdout, "\t%s\n", filename);
    }

    if(readGrid(filename, i * 2, data->computeP1, &data->grid))
      return 1;

    i++;
  }

  /* The grid covers the circumference of the sphere. */
  if(EQUALS(360.0, data->grid.lonStep * data->grid.lonSamples))
    /* Longitude max is infinite. */
    data->grid.lonMax = 1.0e+250;
  else
    /* If not longitude max is that of the last point. */
    data->grid.lonMax = data->grid.lonStep * data->grid.lonSamples;

  data->grid.latMax = data->grid.latStep * data->grid.latSamples;

  return 0;
}











/*
// ///////////////////////////////////////////////////////////////////////////
// Deallocates a grid.
//
// Parameters:
//   grid:   grid to deallocate
*/
static void deleteGrid(gridDsc* grid)
{
  int i, j;

  if(grid->data)
  {
    for (i = 0; i < grid->lonSamples; i++)
    {
      if(grid->data[i])
      {
	for (j = 0; j < grid->latSamples; j++)
	{
	  free(grid->data[i][j]);
	}
	free(grid->data[i]);
      }
    }
  }
  free(grid->data);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Compute the elapsed time since 01/01/1900 0:0:0.0 in julian centuries.
//
// Parameters:
//   date:	Modified Julian day, in seconds.
//
// Return value:
//   The julian centuries.
*/
static double julianCenturies(const double date)
{
  int j;
  int m;
  int y;
  int h;
  int min;
  int sec;
  int mcSecs;
  int bis;
  int day;
  double result;

  decodeDate(date, &j, &m, &y, &h, &min, &sec, &mcSecs);

  bis  = MIN((y % 4) + 1, 2) - 1;

  day = quantiemes[bis][m - 1] + j;

  result  = (y - 1900.0) * 365.0;
  result += day - 1.0;
  result += ((h * 3600.0) + (min * 60.0) + sec + (mcSecs * 1.0E-6)) / 86400.0;

  return (result + nInt((y - 1901.0)/4.0))/ 36525.0;
}


/*
// ///////////////////////////////////////////////////////////////////////////
// Transformation of the description of the waves (amplitude/phase) into
// complex.
//
// Parameters:
//   tide:		The wave describes the tide, if not the radial load.
//   lat:		size of a
//   a:			Array containing the data to be converted.
//
*/
static void ap2complex(const short tide, const int depth, double* a)
{
  int i;
  double amp;
  double pha;

  for(i = 0; i < depth; i+=2)
  {
    amp = a[i    ];
    pha = a[i + 1];

    if(amp != DV && pha != DV)
    {
      a[i    ] = tide?  amp * cos(-pha * rad): amp * cos(pha * rad);
      a[i + 1] = tide? -amp * sin(-pha * rad): amp * sin(pha * rad);
    }
  }
}


/*
// ///////////////////////////////////////////////////////////////////////////
// Perfom bilinear interpolation at point lon, lat from grids.
//
// Parameters:
//   grid:		Grid containing the data->
//   lat:		Latitude
//   lon:		Longitude
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
static int interp(fesData* data, const double lat, const double lon)
{
  int i;
  int rc	= 0;
  int depth	= data->computeP1? 16: 18;
  double re;
  double im;

  if(!IN_GAP(data->westLon, lon, data->eastLon) ||
     !IN_GAP(data->southLat, lat, data->northLat))
  {
    if(getNearestPoints(&data->grid,
			lat,
			lon,
			depth,
			data->sw,
			data->se,
			data->nw,
			data->ne,
			&data->southLat,
			&data->northLat,
			&data->westLon,
			&data->eastLon,
			&data->inGrid))
      goto onError;

      /* Amplitude/phase -> complex */
      ap2complex(data->shortTide, depth, data->sw);
      ap2complex(data->shortTide, depth, data->se);
      ap2complex(data->shortTide, depth, data->nw);
      ap2complex(data->shortTide, depth, data->ne);
  }

  /* La zone recherché n'est pas dans la grille */
  if(!data->inGrid)
    goto noData;

  /* Interpolation des valeurs sur les 18 grilles du fichier NetCdf */
  for(i= 0; i < depth / 2; i++)
  {
    if(data->sw[i << 1] == DV && data->se[i << 1] == DV &&
       data->nw[i << 1] == DV && data->ne[i << 1] == DV)
    {
      goto noData;
    }
    if(bilinearInterp(data->westLon,	/* X1 */
		      data->eastLon,	/* X2 */
		      data->southLat,	/* Y1 */
		      data->northLat,	/* Y2 */
		      data->sw[i << 1],	/* X1, Y1 */
		      data->se[i << 1],	/* X2, Y1 */
		      data->nw[i << 1],	/* X1, Y2 */
		      data->ne[i << 1],	/* X2, Y2 */
		      lon,
		      lat,
		      &re))
      return 1;

    if(re == DV)
      goto noData;

    if(data->sw[(i << 1) + 1] == DV && data->se[(i << 1) + 1] == DV &&
       data->nw[(i << 1) + 1] == DV && data->ne[(i << 1) + 1] == DV)
    {
      goto noData;
    }

    if(bilinearInterp(data->westLon,		/* X1 */
		      data->eastLon,		/* X2 */
		      data->southLat,		/* Y1 */
		      data->northLat,		/* Y2 */
		      data->sw[(i << 1) + 1],	/* X1, Y1 */
		      data->se[(i << 1) + 1],	/* X2, Y1 */
		      data->nw[(i << 1) + 1],	/* X1, Y2 */
		      data->ne[(i << 1) + 1],	/* X2, Y2 */
		      lon,
		      lat,
		      &im))
      return 1;

    if(im == DV)
      goto noData;

    data->re[i] = re;
    data->im[i] = im;
  }

  data->isData = 1;

  /* infer additional constituents by admittance
  DIURNALS (from Richard Ray perth2 program) */

  /* from Q1 and O1 (0-1) */

  /* 2Q1 */
  data->re[16] = 0.263 * data->re[0] - 0.0252 * data->re[1];
  data->im[16] = 0.263 * data->im[0] - 0.0252 * data->im[1];
  /* sigma1 */
  data->re[17] = 0.297 * data->re[0] - 0.0264 * data->re[1];
  data->im[17] = 0.297 * data->im[0] - 0.0264 * data->im[1];
  /* rho1 */
  data->re[18] = 0.164 * data->re[0] + 0.0048 * data->re[1];
  data->im[18] = 0.164 * data->im[0] + 0.0048 * data->im[1];

  /* from O1 and K1  (1-2) */

  /* M11 */
  data->re[19] = 0.0389 * data->re[1] + 0.0282 * data->re[2];
  data->im[19] = 0.0389 * data->im[1] + 0.0282 * data->im[2];
  /* M12 */
  data->re[20] = 0.0140 * data->re[1] + 0.0101 * data->re[2];
  data->im[20] = 0.0140 * data->im[1] + 0.0101 * data->im[2];
  /* chi1 */
  data->re[21] = 0.0064 * data->re[1] + 0.0060 * data->re[2];
  data->im[21] = 0.0064 * data->im[1] + 0.0060 * data->im[2];
  /* pi1 */
  data->re[22] = 0.0030 * data->re[1] + 0.0171 * data->re[2];
  data->im[22] = 0.0030 * data->im[1] + 0.0171 * data->im[2];
  /* phi1 */
  data->re[23] = -0.0015 * data->re[1] + 0.0152 * data->re[2];
  data->im[23] = -0.0015 * data->im[1] + 0.0152 * data->im[2];
  /* theta1 */
  data->re[24] = -0.0065 * data->re[1] + 0.0155 * data->re[2];
  data->im[24] = -0.0065 * data->im[1] + 0.0155 * data->im[2];
  /* J1 */
  data->re[25] = -0.0389 * data->re[1] + 0.0836 * data->re[2];
  data->im[25] = -0.0389 * data->im[1] + 0.0836 * data->im[2];
  /* OO1 */
  data->re[26] = -0.0431 * data->re[1] + 0.0613 * data->re[2];
  data->im[26] = -0.0431 * data->im[1] + 0.0613 * data->im[2];

  /* P1  from Grenoble admittance code */
  if(data->computeP1)
  {
    data->re[8] =  data->aap1 * data->re[0] + data->bbp1 * data->re[1] + data->ccp1 * data->re[2];
    data->im[8] =  data->aap1 * data->im[0] + data->bbp1 * data->im[1] + data->ccp1 * data->im[2];
  }

  /* SEMI-DIURNAL (from Grenoble to take advantage of 2N2) */

  /* from 2N2 -N2 (3-4) */

  /* eps2 */
  data->re[13] = 0.53285 * data->re[3] - 0.03304 * data->re[4];
  data->im[13] = 0.53285 * data->im[3] - 0.03304 * data->im[4];

  /* from M2 - K2 [5-6] */

  /* eta2 */
  data->re[15] = -0.0034925 * data->re[5] + 0.0831707 * data->re[6];
  data->im[15] = -0.0034925 * data->im[5] + 0.0831707 * data->im[6];

  /* from N2 -M2- K2 by spline admittances [see GRL 18[5]:845-848,1991] */

  /* mu2 */
  data->re[10] = data->aamu2 * data->re[6] + data->bbmu2 * data->re[4] + data->ccmu2 * data->re[5];
  data->im[10] = data->aamu2 * data->im[6] + data->bbmu2 * data->im[4] + data->ccmu2 * data->im[5];
  /* nu2 */
  data->re[9] = data->aanu2 * data->re[6] + data->bbnu2 * data->re[4] + data->ccnu2 * data->re[5];
  data->im[9] = data->aanu2 * data->im[6] + data->bbnu2 * data->im[4] + data->ccnu2 * data->im[5];
  /* lda2 */
  data->re[14] = data->aalda2 * data->re[6] + data->bblda2 * data->re[4] + data->cclda2 * data->re[5];
  data->im[14] = data->aalda2 * data->im[6] + data->bblda2 * data->im[4] + data->cclda2 * data->im[5];
  /* L2 */
  data->re[11] = data->aal2 * data->re[6] + data->bbl2 * data->re[4] + data->ccl2 * data->re[5];
  data->im[11] = data->aal2 * data->im[6] + data->bbl2 * data->im[4] + data->ccl2 * data->im[5];
  /* T2 */
  data->re[12] = data->aat2 * data->re[6] + data->bbt2 * data->re[4] + data->cct2 * data->re[5];
  data->im[12] = data->aat2 * data->im[6] + data->bbt2 * data->im[4] + data->cct2 * data->im[5];

  goto onTerminate;

onError:
  rc = 0;

noData:
  data->isData = 0;

onTerminate:
  return rc;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Computes the long-period equilibrium ocean tides.
// Processing logic -
//   Fifteen tidal spectral lines from the Cartwright-Tayler-Edden
//   tables are summed over to compute the long-period tide.
//
// Technical references -
//   Cartwright & Tayler, Geophys. J. R.A.S., 23, 45, 1971.<br>
//   Cartwright & Edden, Geophys. J. R.A.S., 33, 253, 1973.<br>
//
// Parameters:
//   ts:	Julian day, in seconds, denoting time at which tide is
//		to be computed.
//   lat:	Latitude in degrees (positive north) for the position at which
//		tide is computed.
//   tlp:	Computed long-period tide, in centimeters.
*/
void lpeqmt(const double ts, const double lat, double* tlp)
{
  static const double	tc	= 4043174400.0;
  static const double	phc[4]	= { 290.210, 280.120, 274.350, 343.510 };
  static const double	dpd[4]	= { 13.17639650, 0.98564730, 0.11140410, 0.05295390 };
  double		td;
  double		ph;
  double		shpn[4];
  double		zlp;
  double		s;
  double		t	= (ts + 33282.0) * 86400.0;
  int			i;




  /* Compute 4 principal mean longitudes in radians at time TD */
  td = (t - tc) / 86400.0;
  for(i = 0; i < 4; i++)
  {
    ph = phc[i] + (td * dpd[i]);
    shpn[i] = rad * fmod(ph, 360.0);
  }

  /* Assemble long-period tide potential from 15 CTE terms > 1 mm.
  Nodal term is included but not the constant term. */
  zlp = 2.790 * cos(shpn[3])
    -0.490 * cos(shpn[1] - (283.0* rad))
    -3.100 * cos(2.00*shpn[1]);

  ph = shpn[0];
  zlp = zlp - 0.670 * cos(ph - 2.0*shpn[1] + shpn[2])
    - (3.520 - 0.460*cos(shpn[3])) * cos(ph - shpn[2]);
  ph = ph + shpn[0];

  zlp = zlp - 6.660 * cos(ph)
    - 2.760 * cos(ph + shpn[3])
    - 0.260 * cos(ph + 2.0*shpn[3])
    - 0.580 * cos(ph - 2.0*shpn[1])
    - 0.290 * cos(ph - 2.0*shpn[2]);
  ph = ph + shpn[0];
  zlp = zlp - 1.270 * cos(ph - shpn[2])
    - 0.530 * cos(ph - shpn[2] + shpn[3])
    - 0.240 * cos(ph - 2.0 * shpn[1] + shpn[2]);
  /* Multiply by gamma_2 * sqrt(5/4 pi) * P20(lat) */
  s = sin(lat * rad);
  *tlp = 0.4370 * zlp * ((1.50 * SQR(s)) - 0.50);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Init data.
//
//   handle:	  Handle of the internal data.
//   model:	  Tide model (fes99 or fes2002)
//   radialTide:  1 computes radial tide otherwise computes short tide.
//   verbose:	  verbose mode (1/0)
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
int newHandle(void** handle,
	      const tideModel model,
	      const short radialTide,
	      const short verbose)
{
  static double	freq[NB_WAVE] = {
      13.39866087990,	/* Q1		*/
      13.94303558000,	/* O1		*/
      15.04106864000,	/* K1		*/
      27.89535481990,	/* 2N2		*/
      28.43972952010,	/* N2		*/
      28.98410422000,	/* M2		*/
      30.08213728000,	/* K2		*/
      30.00000000000,	/* S2		*/
      14.95893136000,	/* P1		*/
      28.51258314000,	/* Nu2		*/
      27.96820844000,	/* Mu2		*/
      29.52847892000,	/* L2		*/
      29.95893332010,	/* T2		*/
      27.4238337,	/* Eps2		*/
      29.4556253, 	/* Lambda2	*/
      30.6265120, 	/* Eta2		*/
      12.8542862,	/* 2Q1		*/
      12.9271398, 	/* Sigma1	*/
      13.4715145, 	/* Ro1		*/
      14.4966939, 	/* M11		*/
      14.4874103, 	/* M12		*/
      14.5695476, 	/* Ki1		*/
      14.9178647, 	/* Pi1		*/
      15.1232059, 	/* Phi1		*/
      15.5125897, 	/* Teta1	*/
      15.5854433, 	/* J1		*/
      16.1391017, 	/* OO1		*/
  };

  static short code[NB_WAVE] = {
      27,		/* Q1		*/
      1,		/* O1		*/
      3,		/* K1		*/
      5,		/* 2N2		*/
      7,		/* N2		*/
      9,		/* M2		*/
      14,		/* K2		*/
      13,		/* S2		*/
      2,		/* P1		*/
      8,		/* Nu2		*/
      6,		/* Mu2		*/
      11,		/* L2		*/
      12,		/* T2		*/
      60,		/* Eps2		*/
      61,		/* Lambda2	*/
      64,		/* Eta2		*/
      65,		/* 2Q1		*/
      66,		/* Sigma1	*/
      67,		/* Ro1		*/
      68,		/* M11		*/
      69,		/* M12		*/
      70,		/* Ki1		*/
      71,		/* Pi1		*/
      72,		/* Phi1		*/
      73,		/* Teta1	*/
      74,		/* J1		*/
      75,		/* OO1		*/
  };

  int		i;
  double	frbar1;
  double	deno1;
  double	ck;
  double	sk;
  double	cn;
  double	sn;
  double	deno;
  double	cnu2;
  double	snu2;
  double	cmu2;
  double	smu2;
  double	cl2;
  double	sl2;
  double	ct2;
  double	st2;
  double	clda2;
  double	slda2;
  fesData*	data;

  /* Coefficient of the tidal potential (used in spline admittances)
  spline admittances (see GRL 18(5):845-848,1991) */
  #define ALK2		0.1149327
  #define ALN2		0.1758941
  #define ALM2		0.9085024
  #define ALNU2		0.03303
  #define ALMU2		0.02777
  #define ALL2		0.0251
  #define ALLDA2	0.0066
  #define ALT2		0.0247766
  #define ALQ1		0.073017
  #define ALO1		0.3771366
  #define ALK1		0.5300728
  #define ALP1		0.1750754

  /* internal index of the constituen */
  #define IQ1		0
  #define IO1		1
  #define IK1		2
  #define IM2		5
  #define IK2		6
  #define IN2		4
  #define INU2		9
  #define IMU2		10
  #define IL2 		11
  #define IT2 		12
  #define ILDA2		14
  #define IP1		8

  /* Allocate handle */
  if((data = calloc(1, sizeof(fesData))) == NULL)
  {
    fprintf(stderr,
      "Not enought memory to allocate %d bytes.\n",
      sizeof(fesData));
    return 1;
  }

  /* Init tide properties */
  data->computeP1 = model == fes99? 1: 0;
  data->shortTide = !radialTide;
  data->verbose   = verbose;

  if(data->verbose)
  {
    fprintf(stdout, "Initialization of the workspace for computing %s.\n",
    	data->shortTide? "pure oceanic tide": "radial loading tide");
  }

  /* compute the coefficient for P1 as a linear admittance with Q1,O1 K1
  (linear regression) */
  if(data->computeP1)
  {
    frbar1	= 1.0/3.0 * (freq[IQ1] + freq[IO1] + freq[IK1]);
    deno1	= SQR(frbar1) - 1.0 / 3.0 * (SQR(freq[IQ1]) + SQR(freq[IO1]) + SQR(freq[IK1]));

    data->aap1  = ALP1 / 3.0 / ALQ1 *
      (1.0 - (freq[IP1] - frbar1) * (freq[IQ1] - frbar1) / deno1);

    data->bbp1  = ALP1 / 3.0 / ALO1 *
      (1.0 - (freq[IP1] - frbar1) * (freq[IO1] - frbar1) / deno1);

    data->ccp1  = ALP1 / 3.0 / ALK1 *
      (1.0 - (freq[IP1] - frbar1) * (freq[IK1] - frbar1) / deno1);
  }

  ck		= cos(2.0 * pi * 2.0 * (freq[IK2] - freq[IM2]) / 15.0);
  sk		= sin(2.0 * pi * 2.0 * (freq[IK2] - freq[IM2]) / 15.0);

  cn		= cos(2.0 * pi * 2.0 * (freq[IN2] - freq[IM2]) / 15.0);
  sn		= sin(2.0 * pi * 2.0 * (freq[IN2] - freq[IM2]) / 15.0);
  deno		= sk * (cn - 1.0) - sn * (ck - 1.0);

  cnu2		= cos(2.0 * pi * 2.0 * (freq[INU2] - freq[IM2]) / 15.0);
  snu2		= sin(2.0 * pi * 2.0 * (freq[INU2] - freq[IM2]) / 15.0);

  cmu2		= cos(2.0 * pi * 2.0 * (freq[IMU2] - freq[IM2]) / 15.0);
  smu2		= sin(2.0 * pi * 2.0 * (freq[IMU2] - freq[IM2]) / 15.0);

  cl2		= cos(2.0 * pi * 2.0 * (freq[IL2] - freq[IM2]) / 15.0);
  sl2		= sin(2.0 * pi * 2.0 * (freq[IL2] - freq[IM2]) / 15.0);

  ct2		= cos(2.0 * pi * 2.0 * (freq[IT2] - freq[IM2]) / 15.0);
  st2		= sin(2.0 * pi * 2.0 * (freq[IT2] - freq[IM2]) / 15.0);

  clda2		= cos(2.0 * pi * 2.0 * (freq[ILDA2] - freq[IM2]) / 15.0);
  slda2		= sin(2.0 * pi * 2.0 * (freq[ILDA2] - freq[IM2]) / 15.0);

  data->aamu2	=
    (-sn * cmu2 + (cn - 1.0) * smu2 + sn) / deno/ ALK2 * ALMU2;
  data->aanu2	=
    (-sn * cnu2 + (cn - 1.0) * snu2 + sn) / deno / ALK2 * ALNU2;
  data->aal2	=
    (-sn * cl2  + (cn - 1.0) * sl2  + sn) / deno / ALK2 * ALL2;
  data->aat2	=
    (-sn * ct2  + (cn - 1.0) * st2  + sn) / deno / ALK2 * ALT2;
  data->aalda2  =
    (-sn * clda2+ (cn - 1.0) * slda2+ sn) / deno / ALK2 * ALLDA2;

  data->bbmu2	=
    (sk * cmu2 - (ck - 1.0) * smu2 - sk) / deno / ALN2 * ALMU2;
  data->bbnu2	=
    (sk * cnu2 - (ck - 1.0) * snu2 - sk) / deno / ALN2 * ALNU2;
  data->bbl2	=
    (sk * cl2  - (ck - 1.0) * sl2  - sk) / deno / ALN2 * ALL2;
  data->bbt2	=
    (sk * ct2  - (ck - 1.0) * st2  - sk) / deno / ALN2 * ALT2;
  data->bblda2  =
    (sk * clda2- (ck - 1.0) * slda2- sk) / deno/ ALN2 * ALLDA2;

  data->ccmu2	=
    (-(sk - sn) * cmu2 + (ck - cn) * smu2 + sk * cn - sn * ck)
    / deno / ALM2 * ALMU2;
  data->ccnu2	=
    (-(sk - sn) * cnu2 + (ck - cn) * snu2 + sk * cn - sn * ck)
    / deno / ALM2 * ALNU2;
  data->ccl2	=
    (-(sk - sn) * cl2  + (ck - cn) * sl2  + sk * cn - sn * ck)
    / deno / ALM2 * ALL2;
  data->cct2	=
    (-(sk - sn) * ct2  + (ck - cn) * st2  + sk * cn - sn * ck)
    / deno / ALM2 * ALT2;
  data->cclda2  =
    (-(sk - sn) * clda2+ (ck - cn) * slda2+ sk * cn - sn * ck)
    / deno / ALM2 * ALLDA2;

  for(i = 0; i < NB_WAVE; i++)
  {
    data->freq[i]	= freq[i];
    data->num[i]	= code[i];
    data->freq[i]	= data->freq[i] * rad;
  }

  initCorrections(data, 0.0);

  *handle = data;

  return 0;
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
int fesCore(void* handle,
	    const double lat,
	    const double lon,
	    const double time,
	    double* tide)
{
  int	    i;
  double    phi;
  double    t1	    = julianCenturies(time * 86400.0);
  double    delta;
  fesData*  data    = (fesData*)handle;

  delta = (t1 - data->tNodal) * 3.6525E+04 * 24.0;

  if(fabs(delta) > 24.0)
  {
    initCorrections(handle, t1);
    delta = 0.0;
  }

  if(interp(handle, lat, lon))
    return 1;

  if(data->isData)
  {
    *tide = 0.0;

    for(i = 0; i < NB_WAVE; i++)
    {
      phi = fmod(data->freq[i] * delta + data->v0u[i], 2.0 * pi);

      if (phi < 0.0)
	phi = phi + 2.0 * pi;

      *tide += data->f[i] * (data->re[i] * cos(phi) + data->im[i] * sin(phi));
    }
  }
  else
  {
    *tide = DV;
  }

  return 0;
}

/*
// ///////////////////////////////////////////////////////////////////////////
// Deallocates the internal data.
//
// Parameters:
//   handle:	Handle of the internal data.
*/
void freeHandle(void* handle)
{
  fesData* data = (fesData*) handle;

  if(data->verbose)
  {
    fprintf(stdout, "Release workspace for computing %s.\n",
    	data->shortTide? "pure oceanic tide": "radial loading tide");
  }

  deleteGrid(&data->grid);
  free(data);
}
/* eof fescore.c */
