#include "fescore.h"
#ifndef FESTIDE_H
#define FESTIDE_H

/*
// ///////////////////////////////////////////////////////////////////////////
// TODO
//
// Parameters:
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
int initFes(const char* const dir, const tideModel model)
;

/*
// ///////////////////////////////////////////////////////////////////////////
// TODO
//
// Parameters:
//
// Return value:
//   1 is returned to indicate an error otherwise 0.
*/
void freeFes(void)
;


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
;
#endif
/* eof festide.h */
