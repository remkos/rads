/* ######################################################################
 *  Main routines for the FES prediction software.
 * 
 *  File      : fes.c
 *  Developer : CLS
 *  Version   : 1.1
 *  Date      : 6 oct 2004
 *  
 * ######################################################################
 */

#include "fes-int.h"

#ifdef BINARY
#ifdef WIN32
/* NTFS file system */
static const char* tideFmt  = "%s\\%s_fes%s.bin";
static const char* loadFmt  = "%s\\%s_drfes%s.bin";
#else
/* UNIX file system */
static const char* tideFmt  = "%s/%s_fes%s.bin";
static const char* loadFmt  = "%s/%s_drfes%s.bin";
#endif
#else
#ifdef WIN32
/* NTFS file system */
static const char* tideFmt  = "%s\\%s_fes%s.asc";
static const char* loadFmt  = "%s\\%s_drfes%s.asc";
#else
/* UNIX file system */
static const char* tideFmt  = "%s/%s_fes%s.asc";
static const char* loadFmt  = "%s/%s_drfes%s.asc";
#endif
#endif

/* constants */
static const char* waveName[] =
{
    "Q1", "O1", "K1", "2N2", "N2", "M2", "K2", "S2", "P1", NULL
};

static const char* modelName[] =
{
    "99", "2002", "2004", NULL
};





/*
// ///////////////////////////////////////////////////////////////////////////
// initializes the computation of the tide.
//
//   fesData:	Internal data handle, which is a pointer to a fesData
//		structure containing information about the tide computation.
//   mode:	Computation mode. If mode is equals to FES_SHORT_TIDE,
//		fesCore computes the short tide otherwise she computes the
//		radial tide.
//   dir:    	Directory containing grids
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.
*/
int fesNew(fesData** fes,
	   const int shortTide,
	   const tideModel model,
	   const char* const dir)
{
    char   filename[MAX_PATH];
    int    ix;
    int    rc;
#ifndef IO
    FILE*  handle;
#endif

    assert(*fes == NULL);

    /* Allocate handle */
    if ( (*fes = calloc(1, sizeof(fesData))) == NULL )
	return FES_NO_MEMORY;

    (*fes)->grid.depth     = model == fes99? DEPTH_FES99: DEPTH_FES2002;
    (*fes)->grid.shortTide = shortTide;
    (*fes)->computeP1	   = model == fes99;
    (*fes)->tNodal	   = 0;

#ifdef IO
    /* Allocate handle */
    if ( ((*fes)->grid.handle = calloc((*fes)->grid.depth,
	sizeof(FILE*))) == NULL )
	return FES_NO_MEMORY;
#else
    /* Allocate buffer */
    if ( ((*fes)->grid.buffer = calloc((*fes)->grid.depth,
	sizeof(fComplex*))) == NULL )
	return FES_NO_MEMORY;
#endif


    for ( ix = 0; ix < (*fes)->grid.depth; ix++ )
    {
        /*printf ("Reading grid for: %s\n", waveName[ix]);*/
	/* ignore P1 for fes99 */
	if( strcmp(waveName[ix], "P1") || model != fes99 )
	{
	    sprintf(filename,
		(*fes)->grid.shortTide ? tideFmt: loadFmt,
		dir,
		waveName[ix],
                modelName[model]);
#ifdef IO
	    if ( ((*fes)->grid.handle[ix] = fopen(filename, "rb")) == NULL )
		return FES_ACCESS_ERROR;

	    rc = readHeader((*fes)->grid.handle[ix], ix, &(*fes)->grid);
#else
	    if ( (handle = fopen(filename, "rb")) == NULL )
		return FES_ACCESS_ERROR;

	    rc = readHeader(handle, ix, &(*fes)->grid);
	    fclose(handle);
#endif
	    if ( rc != FES_SUCCESS )
		return rc;

	}
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
#ifdef IO
	    if ( fes->grid.handle[ix] )
	        fclose(fes->grid.handle[ix]);
        free(fes->grid.handle);
#else
	    free(fes->grid.buffer[ix]);
        free(fes->grid.buffer);
#endif
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

    admittance(fes->computeP1, fes->waves);

    if ( fes->isData )
    {
	*h = 0.0;

	for ( i = 0; i < N_WAVES; i++ )
	{
	    phi = fmod( fes->waves[i].freq * delta + fes->waves[i].v0u, 2.0 * M_PI );

	    if ( phi < 0.0 )
		phi = phi + 2.0 * M_PI;

	    *h += fes->waves[i].f *
		(fes->waves[i].cplx.re * cos(phi) + fes->waves[i].cplx.im * sin(phi));
	}

	if ( fes->grid.shortTide )
	{
	    lpeqmt(time, lat, hLp);
	}
	else
	    *hLp = 0;
    }
    else
	return FES_NO_DATA;

    return FES_SUCCESS;
}
