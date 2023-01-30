#ifndef FESCORE_H
#define FESCORE_H


/*
// ///////////////////////////////////////////////////////////////////////////
// Definitions of the models of known tides.
*/
typedef enum {fes99 = 0, fes2002} tideModel;




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
;






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
;





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
;





/*
// ///////////////////////////////////////////////////////////////////////////
// Driver for tidal computation.
//
// Parameters:
//   handle:	        handle of the internal data.
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
;


/*
// ///////////////////////////////////////////////////////////////////////////
// Deallocates the internal data.
//
// Parameters:
//   handle:	Handle of the internal data.
*/
void freeHandle(void* handle)
;

#endif
/* eof fescore.h */
