/*
 *  Main routines for the FES prediction software.
 * 
 *  File      : fes.c
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.4
 *  Date      : 11 May 2005
 *  
 *  1.2: add M4, S1 (Ray), Mf, Mm, Mtm, MSqm
 *  1.3: remove conditional compilation
 *  1.4: correct computation of LP
 */

#include "fes-int.h"

#ifdef WIN32
/* NTFS file system */
static const char* tideFmt  = "%s\\%s_fes%s.%s";
static const char* loadFmt  = "%s\\%s_drfes%s.%s";
#else
/* UNIX file system */
static const char* tideFmt  = "%s/%s_fes%s.%s";
static const char* loadFmt  = "%s/%s_drfes%s.%s";
#endif

/* constants */
static const char* waveName[] =
{
    "Q1", "O1", "K1", "2N2", "N2", "M2", "K2", "S2", "P1", "M4", "S1", "Mf", "Mm", "Mtm", "MSqm", NULL
};

static const char* modelName = "2004";




/*
// ///////////////////////////////////////////////////////////////////////////
// initializes the computation of the tide.
//
//   fes:	Internal data handle, which is a pointer to a fesData
//		structure containing information about the tide computation.
//   tide: 	Computation mode. If mode is equals to FES_TIDE, 
//		fesCore computes the tide otherwise she computes the radial
//		tide.
//   mode:	One  of  FES_MEM, FES_IO, FES_ASCII which request loading 
//		grids into memory, direct access from binary grids or
//		loading ASCII grids into memory
//   dir:    	Directory containing grids
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.
*/
int fesNew(fesData** fes,
	   const int tide,
	   const int mode,
	   const char* const dir)
{
    char   filename[MAX_PATH];
    int    ix;
    int    rc;
    FILE*  handle		= NULL;

    assert(*fes == NULL);

    /* Allocate handle */
    if ( (*fes = calloc(1, sizeof(fesData))) == NULL )
	return FES_NO_MEMORY;

    (*fes)->grid.depth	= tide? DEPTH_TIDE: DEPTH_RADIAL;
    (*fes)->grid.tide	= tide;
    (*fes)->tNodal	= 0;

    switch(mode)
    {
    /* Allocate handle */
    case FES_IO:
	if ( ((*fes)->grid.handle = calloc((*fes)->grid.depth,
	    sizeof(FILE*))) == NULL )
	    return FES_NO_MEMORY;
        break;
    /* Allocate buffer */
    case FES_MEM:
    case FES_ASCII:
	if ( ((*fes)->grid.buffer = calloc((*fes)->grid.depth,
	    sizeof(fComplex*))) == NULL )
	    return FES_NO_MEMORY;
        break;
    default:
      return FES_MODE_ERROR;
    }

    for ( ix = 0; ix < (*fes)->grid.depth; ix++ )
    {
	sprintf(filename,
	    (*fes)->grid.tide ? tideFmt: loadFmt,
	    dir,
	    waveName[ix],
            modelName,
	    mode == FES_ASCII? "asc": "bin");
        /* direct access */
	if( (*fes)->grid.handle != NULL )
	{
	    if ( ((*fes)->grid.handle[ix] = fopen(filename, "rb")) == NULL )
		return FES_ACCESS_ERROR;

	    rc = readHeader((*fes)->grid.handle[ix], mode, ix, &(*fes)->grid);
	}
	/* loading grid into memory */
	else
	{
	    if ( (handle = fopen(filename, "rb")) == NULL )
		return FES_ACCESS_ERROR;

	    rc = readHeader(handle, mode, ix, &(*fes)->grid);
	    fclose(handle);
	}
	if ( rc != FES_SUCCESS )
	    return rc;
    }

    /* The grid covers the circumference of the sphere. */
    if ( EQUALS(360.0, (*fes)->grid.lonStep * (*fes)->grid.lonSamples) )
	/* LongitudLe max is infinite. */
	(*fes)->grid.lonMax = 1.0e+250;
    else
	/* If not longitude max is that of the last point. */
	(*fes)->grid.lonMax = getValue((*fes)->grid.lonSamples - 1,
    				   (*fes)->grid.lonMin,
    				   (*fes)->grid.lonStep);

    (*fes)->grid.latMax = getValue((*fes)->grid.latSamples - 1,
    				   (*fes)->grid.latMin,
    				   (*fes)->grid.latStep);

    initAdmittance((*fes)->waves);
    initCorrections(0,(*fes)->waves);

    return FES_SUCCESS;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Frees a fesData structure from memory.
//
// Parameters:
//   fes:	Pointer to the fesData structure that you want to free.
*/
void fesDelete(fesData* fes)
{
    int     ix;
    
    if( fes != NULL )
    {
        for ( ix = 0; ix < fes->grid.depth; ix++ )
	{
	    if( fes->grid.handle != NULL )
		if ( fes->grid.handle[ix] != NULL )
	            fclose(fes->grid.handle[ix]);

	    if( fes->grid.buffer != NULL )
		free(fes->grid.buffer[ix]);
	    
	}
        free(fes->grid.handle);
        free(fes->grid.buffer);
        free(fes);
    }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Driver for tidal computation.
//
// Parameters:
//   fes:	        Internal data handle, which is a pointer to a fesData
//			structure containing information about the tide
//			computation.
//   lat		Latitude in degrees (positive north) for the position
//			at which tide is computed.
//   lon		Longitude in degrees (positive north) for the position
//			at which tide is computed.
//   time		Julian day (days since 1950-01-01 00:00:00.000 UTC).
//   h			Computed height of the diurnal and semi-diunral
//                      constituents of the tidal spectrum (in centimeters)
//   hLp		Computed height of the long period wave
//                      constituents of the tidal spectrum (in centimeters)
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.
*/
int fesCore(fesData* fes,
	    const double lat,
	    const double lon,
	    const double time,
	    double* h,
	    double* hLp)
{
    int       rc;
    int       i;
    double    phi;
    double    t1      = julianCenturies(time);
    double    delta;

    delta = (t1 - fes->tNodal) * 3.6525E+04 * 24.0;

    if ( fabs(delta) > 24.0 )
    {
	initCorrections(t1, fes->waves);
	fes->tNodal = t1;
	delta   = 0.0;
    }

    rc = interp(fes, lat, lon);
    if ( rc != FES_SUCCESS )
	return rc;

    admittance(fes->waves);

    if ( fes->isData )
    {
	*h = 0.0;

	if ( fes->grid.tide )
	{
	    lpeMinus4Waves(time, lat, hLp);
	}
	else
	    *hLp = 0;

	for ( i = 0; i < N_WAVES; i++ )
	{
	    double	tide;

            phi = fmod( fes->waves[i].freq * delta + fes->waves[i].v0u, 2.0 * M_PI );

	    if ( phi < 0.0 )
		phi = phi + 2.0 * M_PI;
                
            tide = fes->waves[i].f *
		  (fes->waves[i].cplx.re * cos(phi) + fes->waves[i].cplx.im * sin(phi));

	    if( fes->waves[i].type == SP_TIDE )
              *h   += tide;
            else
              *hLp += tide;
	}
    }
    else
	return FES_NO_DATA;

    return FES_SUCCESS;
}
