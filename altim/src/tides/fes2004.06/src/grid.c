/*
 *  Management of grids.
 * 
 *  File      : grid.c
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.3
 *  Date      : 13 January 2005
 *  
 *  1.2: add M4, S1 (Ray), Mf, Mm, Mtm, MSqm
 *  1.3: remove conditional compilation
 */

#include "fes-int.h"

#define N                   30

/*
// ///////////////////////////////////////////////////////////////////////////
// This function swaps a variable
//
// Parameters:
//   a:		pointer to the variable to swap
//   b:		pointer to the variable to swap
*/
static void swap(void** a, void** b)
{
    register void *tmp = *a;

    *a = *b;
    *b = tmp;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// function swaps a "size" byte variable, i.e. it converts their
// representation between little and big endian (in either direction).
//
// Parameters:
//   p:		pointer to the variable to swap
//   n:		size of the variable
*/
static void swapByte(void *p, int n)
{
    register char* begin        = (char*)p;
    register char* end          = ((char*)p) +(n - 1);
    register char  temp;

    while ( begin < end )
    {
	temp = *begin;
	*(begin++) = *end;
	*(end--)   = temp;
    }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Verifie the system integer variables notation.
//
// Return value:
//   Returns 1 for little endian else 0 for big endian notation.
*/
static int isLittleEndian(void)
{
    static short a = 1;
    static char *p = (char *)&a;

    return(*p ? 1 : 0);
}




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

    while ( result >= ((base + 360.0) - EPSILON) )
	result -= 360.0;
    while ( result < base - EPSILON )
	result += 360.0;
    if ( fabs(result - base) <= EPSILON )
	result = base;

    return result;
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
    return nInt((value - min) / step);
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
double getValue(const int idx, const double min, const double step)
{
    return min + (step * idx);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Reading a value in the bmg file.
//
// Parameters:
//   grid:	Grid properties
//   iLon:	Index of the longitude in the grid.
//   iLat:	Index of the latitude in the grid.
//   iDepth:	Index of the grid file.
//   value:	value readed.
//
// Return value:
//   If an error occurs, the return value is FES_INPUT_ERROR, otherwise
//   FES_SUCCESS.
*/
static int readGridValue(gridDsc* grid,
			 const int iLon,
			 const int iLat,
			 const int iDepth,
			 dComplex* value)
{
    fComplex    z;

    /* reading value from grid */
    if( grid->handle != NULL )
    {
	long    offset;

	offset = 28 + (iLat * grid->lonSamples + iLon) * sizeof(fComplex);

	/* Sets the file position to read the value */
	if ( fseek(grid->handle[iDepth], offset, SEEK_SET) ==  -1L )
	    return FES_INPUT_ERROR;

	if ( fread(&z, sizeof(fComplex), 1, grid->handle[iDepth]) != 1 )
	    return FES_INPUT_ERROR;

	/* Swap bytes ? */
	if ( !grid->littleEndian )
	{
	    swapByte(&z, sizeof(fComplex));
	    swap((void**)&z.re, (void**)&z.im);
	}
    }
    /* reading values from memory */
    else
    {
	z   = grid->buffer[iDepth][iLat * grid->lonSamples + iLon];
    }

    /* set result to default value ? */
    if ( z.re == grid->undef || z.im == grid->undef )
    {
	value->re = DV;
	value->im = DV;
    }
    else
    {
	value->re = z.re;
	value->im = z.im;
    }
    return FES_SUCCESS;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Reading the points closest to a position.
//
// Parameters:
//   grid:	Grid containing the data
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
//   If an error occurs, the return value is FES_INPUT_ERROR, otherwise
//   FES_SUCCESS.
*/
static int getNearestPoints(gridDsc* grid,
			    const double lat,
			    const double lon,
			    const int depth,
			    dComplex* sw,
			    dComplex* se,
			    dComplex* nw,
			    dComplex* ne,
			    double* southLat,
			    double* northLat,
			    double* westLon,
			    double* eastLon,
			    int* inGrid)
{
    int     rc;
    double  lonNorm = normalizeLongitude(grid->lonMin, lon);

    assert( depth == grid->depth );

    /* Check if asked position is in the grid */
    if ( (!IN_GAP(grid->latMin, lat, grid->latMax)) ||
	 (!IN_GAP(grid->lonMin, lonNorm, grid->lonMax)) )
    {
	*inGrid = 0;
    }
    else
    {
	int iLat1   = getIndex(grid->latMin, lat, grid->latStep);
	int iLon1   = getIndex(grid->lonMin, lonNorm, grid->lonStep);
	int iLat2;
	int iLon2;
	int iDepth;

	if ( lat >= grid->latMax )
	{
	    iLat2 = iLat1;
	    iLat1--;
	}
	else
	    iLat2 = iLat1 + 1;

	if ( lonNorm >= grid->lonMax )
	{
	    iLon2 = iLon1;
	    iLon1--;
	}
	else
	    iLon2 = iLon1 + 1;

	*inGrid = 1;

	*southLat = getValue(iLat1, grid->latMin, grid->latStep);
	*northLat = getValue(iLat2, grid->latMin, grid->latStep);
	*westLon  = getValue(iLon1, grid->lonMin, grid->lonStep);
	*eastLon  = getValue(iLon2, grid->lonMin, grid->lonStep);

	if ( *westLon != lon )
	{
	    double gap = *westLon - *eastLon;

	    lonNorm  = normalizeLongitude(lon, *eastLon);
	    *eastLon = lonNorm;
	    *westLon = lonNorm + gap;
	}

	iLon1 %= grid->lonSamples;
	iLon2 %= grid->lonSamples;

	for ( iDepth = 0; iDepth < depth; iDepth++ )
	{
	    rc = readGridValue(grid, iLon1, iLat1, iDepth, &sw[iDepth]);
	    if ( rc != FES_SUCCESS )
		return rc;

	    rc = readGridValue(grid, iLon2, iLat1, iDepth, &se[iDepth]);
	    if ( rc != FES_SUCCESS )
		return rc;

	    rc = readGridValue(grid, iLon1, iLat2, iDepth, &nw[iDepth]);
	    if ( rc != FES_SUCCESS )
		return rc;

	    rc = readGridValue(grid, iLon2, iLat2, iDepth, &ne[iDepth]);
	    if ( rc != FES_SUCCESS )
		return rc;
	}
    }
    return FES_SUCCESS;
}




/*
// ///////////////////////////////////////////////////////////////////////////
// Read formatted data from a line.
//
// Parameters:
//   line:	line to read
//   sizeOf	size of values
//   values:	The data read.
//
// Return value:
//   returns the number of input fields successfully scanned, converted
//   and stored
*/
static int readValues(char* line,
		      const int sizeOf,
		      float* values)
{
    int   ix      = 0;
    char* token   = strtok(line, " \t");

    while(token)
    {
	if( ix > sizeOf )
	    return -1;

	if(sscanf(token, "%f", &values[ix]) != 1)
	    return ix;
	ix += 1;

	token  = strtok(NULL, " \t");
    }

    return ix;
}




/*
// ///////////////////////////////////////////////////////////////////////////
// Reads header file.
//
// Parameters:
//   handle	Handle to the current grid.
//   mode:	One  of  FES_MEM, FES_IO, FES_ASCII which request loading 
//		grids into memory, direct access from binary grids or
//		loading ASCII grids into memory
//   iDepth	current depth
//   grid:	Grid containing the data readed.
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.
*/
int readHeader(FILE* handle,
	       const int mode,
	       const int iDepth,
	       gridDsc* grid)
{
    int     rc;
    int     nX;
    int     nY;
    float   latMin;
    float   lonMin;
    float   latStep;
    float   lonStep;
    float   undef;
    char    buffer[MAX_PATH];
    float   unused;
    float   amp[N];
    float   pha[N];
    int     x;
    int     y;

    /* reading ascii grids */
    if( mode == FES_ASCII )
    {
	/* Reading lat min, max */
	if(fscanf(handle, "%f %f", &lonMin, &unused) != 2)
	    return FES_INPUT_ERROR;

	/* Reading lon min, max */
	if(fscanf(handle, "%f %f", &latMin, &unused) != 2)
	    return FES_INPUT_ERROR;

	/* Reading lon, lat step */
	if(fscanf(handle, "%f %f", &lonStep, &latStep) != 2)
	    return FES_INPUT_ERROR;

	/* Reading nX, nY */
	if(fscanf(handle, "%d %d", &nX, &nY) != 2)
	    return FES_INPUT_ERROR;

	/* The grids should not overlap. */
	nX--;

	/* Reading mask */
	if(fscanf(handle, "%f %f", &undef, &unused) != 2)
	    return FES_INPUT_ERROR;
	/* flush line before read values */
	if(fgets(buffer, MAX_PATH, handle) == NULL)
	    return FES_INPUT_ERROR;
    }
    /* reading binary grids */
    else
    {
	grid->littleEndian = isLittleEndian();

	/* Reading lat, lon samples */
	rc = fread(&nX, sizeof(int), 1, handle);
	if ( rc != 1 )
	    return FES_INPUT_ERROR;

	rc = fread(&nY, sizeof(int), 1, handle);
	if ( rc != 1 )
	    return FES_INPUT_ERROR;

	if ( !grid->littleEndian )
	{
	    swapByte(&nX, sizeof(int));
	    swapByte(&nY, sizeof(int));
	}

	/* Reading lat, lon min */
	rc = fread(&lonMin, sizeof(float), 1, handle);
	if ( rc != 1 )
	    return FES_INPUT_ERROR;

	rc = fread(&latMin, sizeof(float), 1, handle);
	if ( rc != 1 )
	    return FES_INPUT_ERROR;

	/* Reading lat, lon step */
	rc = fread(&lonStep, sizeof(float), 1, handle);
	if ( rc != 1 )
	    return FES_INPUT_ERROR;

	rc = fread(&latStep, sizeof(float), 1, handle);
	if ( rc != 1 )
	    return FES_INPUT_ERROR;

	/* Reading default value */
	rc = fread(&undef, sizeof(float), 1, handle);
	if ( rc != 1 )
	    return FES_INPUT_ERROR;

	if ( !grid->littleEndian )
	{
	    swapByte(&lonMin, sizeof(float));
	    swapByte(&latMin, sizeof(float));
	    swapByte(&lonStep, sizeof(float));
	    swapByte(&latStep, sizeof(float));
	    swapByte(&undef, sizeof(float));
	}
    }

    /* First call */
    if ( iDepth == 0 )
    {
	grid->latMin        = latMin;
	grid->lonMin        = lonMin;
	grid->latStep       = latStep;
	grid->lonStep       = lonStep;
	grid->latSamples    = nY;
	grid->lonSamples    = nX;
	grid->undef         = undef;
    }
    /* Check grid definition */
    else if ( grid->latMin     != latMin   || grid->lonMin     != lonMin  ||
	      grid->latStep    != latStep  || grid->lonStep    != lonStep ||
	      grid->latSamples != nY       || grid->lonSamples != nX      ||
	      grid->undef      != undef )
	return FES_GRIDS_CONSTANT;

    /* loading grids into memory */
    if( grid->buffer != NULL )
    {
	int     ix;
	int     size	= nX * nY;

	/* Allocate the current grid */
	if ( (grid->buffer[iDepth] = calloc(size, sizeof(fComplex))) == NULL )
	    return FES_NO_MEMORY;

	/* loading ascii grid */
	if( mode == FES_ASCII )
	{
	    y = 0;

	    /* For all the latitudes. */
	    while( y < nY )
	    {
		x = 0;

		/* For all the longitudes. */
		while( x < nX )
		{
		    /* Reading the 30 values of the amplitude. */
		    if( fgets(buffer, MAX_PATH, handle) == NULL )
			return FES_INPUT_ERROR;
		    if( readValues(buffer, N, amp) != N )
			return FES_INPUT_ERROR;

		    /* Reading the 30 values of the phase. */
		    if( fgets(buffer, MAX_PATH, handle) == NULL )
			return FES_INPUT_ERROR;
		    if( readValues(buffer, N, pha) != N )
			return FES_INPUT_ERROR;

		    for( ix = 0; ix < N; ix++ )
		    {
			if( amp[ix] != undef && pha[ix] != undef )
			{
			    grid->buffer[iDepth][y * nX + x + ix].re = grid->tide?
				 amp[ix] * cos(-pha[ix] * RAD):
				 amp[ix] * cos( pha[ix] * RAD);
			    grid->buffer[iDepth][y * nX + x + ix].im = grid->tide?
				-amp[ix] * sin(-pha[ix] * RAD):
				 amp[ix] * sin( pha[ix] * RAD);
			}
			else
			{
			    grid->buffer[iDepth][y * nX + x + ix].re = undef;
			    grid->buffer[iDepth][y * nX + x + ix].im = undef;
			}
		    }
		    x += N;
		}
		/* Reading of the redundant value for the amplitude. */
		if(fgets(buffer, MAX_PATH, handle) == NULL)
		    return FES_INPUT_ERROR;

		/* Reading of the redundant value for the phase. */
		if(fgets(buffer, MAX_PATH, handle) == NULL)
		    return FES_INPUT_ERROR;
		y++;
	    }
        }
	/* loading binary grid */
	else
	{
	    /* reading all values */
	    rc = fread(grid->buffer[iDepth], sizeof(fComplex), size, handle);
	    if ( rc != size )
		return FES_INPUT_ERROR;

	    /* Swap bytes ? */
	    if ( !grid->littleEndian )
	    {
		for ( ix = 0; ix < size; ix++ )
		{
		    swapByte(&grid->buffer[iDepth][ix], sizeof(fComplex));
		    swap((void*)&grid->buffer[iDepth][ix].re,
			 (void*)&grid->buffer[iDepth][ix].im);
		}
	    }
	}
    }
    return FES_SUCCESS;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Perfom bilinear interpolation at point lon, lat from grids.
//
// Parameters:
//   grid:		Grid containing the data
//   lat:		Latitude
//   lon:		Longitude
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.
*/
int interp(fesData* fes, const double lat, const double lon)
{
    int     rc;
    int     ix;
    dComplex    cplx;

    if ( !IN_GAP(fes->westLon, lon, fes->eastLon) ||
	 !IN_GAP(fes->southLat, lat, fes->northLat) )
    {
	rc = getNearestPoints(&fes->grid,
			      lat,
			      lon,
			      fes->grid.depth,
			      fes->sw,
			      fes->se,
			      fes->nw,
			      fes->ne,
			      &fes->southLat,
			      &fes->northLat,
			      &fes->westLon,
			      &fes->eastLon,
			      &fes->inGrid);
	if ( rc != FES_SUCCESS )
	    return rc;
    }

    /* The zone required is not in the grid */
    if ( !fes->inGrid )
	goto noData;

    /* Interpolation */
    for ( ix= 0; ix < fes->grid.depth; ix++ )
    {
	if ( fes->sw[ix].re == DV && fes->se[ix].re == DV &&
 	     fes->nw[ix].re == DV && fes->ne[ix].re == DV )
	    goto noData;

	bilinearInterp(fes->westLon,    /* X1 */
		       fes->eastLon,    /* X2 */
		       fes->southLat,   /* Y1 */
		       fes->northLat,   /* Y2 */
		       fes->sw[ix].re,  /* X1, Y1 */
		       fes->se[ix].re,  /* X2, Y1 */
		       fes->nw[ix].re,  /* X1, Y2 */
		       fes->ne[ix].re,  /* X2, Y2 */
		       lon,
		       lat,
		       &cplx.re);

	if ( cplx.re == DV )
	    goto noData;

	bilinearInterp(fes->westLon,    /* X1 */
		       fes->eastLon,    /* X2 */
		       fes->southLat,   /* Y1 */
		       fes->northLat,   /* Y2 */
		       fes->sw[ix].im,  /* X1, Y1 */
		       fes->se[ix].im,  /* X2, Y1 */
		       fes->nw[ix].im,  /* X1, Y2 */
		       fes->ne[ix].im,  /* X2, Y2 */
		       lon,
		       lat,
		       &cplx.im);

	if ( cplx.im == DV )
	    goto noData;

	fes->waves[ix].cplx = cplx;
    }

    fes->isData = 1;

    goto onTerminate;

noData:
    fes->isData = 0;

onTerminate:
    return FES_SUCCESS;
}
