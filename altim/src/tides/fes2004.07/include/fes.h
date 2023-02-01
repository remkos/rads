/*
 *  Include file for the FES prediction software.
 * 
 *  File      : fes.h
 *  Developer : CLS
 *  Version   : 1.3
 *  Date      : 13 January 2005
 *  
 *  1.2: add M4, S1 (Ray), Mf, Mm, Mtm, MSqm
 *  1.3: remove conditional compilation
 *  
 *  This software have been tested on Linux platform
 *  It is provided without any guarantees.
 */

#ifndef _FES_H
#define _FES_H

#ifdef __cplusplus
extern "C" {
#endif

/* possible error codes we can be returned */
#define FES_SUCCESS		0x00	/* 0 */
#define FES_INPUT_ERROR		0x01	/* 1 */
#define FES_ACCESS_ERROR        0x02	/* 2 */
#define FES_GRIDS_CONSTANT      0x03	/* 3 */
#define FES_INV_ERROR      	0x04	/* 4 */
#define FES_NO_MEMORY		0x05	/* 5 */
#define FES_NO_DATA		0x06	/* 6 */
#define FES_MODE_ERROR		0x07	/* 7 */

/* possible tide */
#define FES_RADIAL		0
#define FES_TIDE		1
/* possible access */
#define FES_MEM			0
#define FES_IO			1
#define FES_ASCII		2

    /* opaque connection handle */
    typedef struct _fesData fesData;





    /*
    // ///////////////////////////////////////////////////////////////////////////
    // initializes the computation of the tide.
    //
    //   fesData:   Internal data handle, which is a pointer to a fesData 
    //		    structure containing information about the tide computation.
    //   shortTide: Computation mode. If mode is equals to FES_SHORT_TIDE, 
    //		    fesCore computes the short tide otherwise she computes the 
    //		    radial tide.
    //   mode:	    One  of  FES_MEM, FES_IO which request loading grids into
    //		    memory or direct access
    //   dir:	    Directory containing grids
    //
    // Return value:
    //   Returns FES_SUCCESS if the operation completed successfully or an error
    //   code if a problem occurred.
    */
    int fesNew(fesData** fes,
	       const int shortTide,
	       const int mode,
	       const char* const dir);





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
    //   h		Computed height of the diurnal and semi-diunral
    //                  constituents of the tidal spectrum (in centimeters)
    //   hLp		Computed height of the long period wave
    //                  constituents of the tidal spectrum (in centimeters)
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
		double* hLp);





    /*
    // ///////////////////////////////////////////////////////////////////////////
    // Frees a fesData structure from memory.
    //
    // Parameters:
    //   fes:	Pointer to the fesData structure that you want to free.
    */
    void fesDelete(fesData* fes);





    /*
    // ///////////////////////////////////////////////////////////////////////////
    // The function returns a pointer to the error message.
    //
    // Parameters:
    //   err:	Error code that you want interpreted into an error message.
    //
    // Return value:
    //   A pointer to the error message.
    */
    char* fesErr2String( const int err );

#ifdef __cplusplus
}
#endif
#endif /* _FES_H */
