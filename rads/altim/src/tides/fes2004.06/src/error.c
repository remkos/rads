/*
 *  Management of errors.
 * 
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.2
 *  Date      : 13 January 2005
 *  
 *  1.2: remove conditional compilation
 *  
 */

#include "fes-int.h"

/*
* This structure represents FES errors messages.
*/
struct _fesErr
{
    int   code;
    char* reason;
};

static struct _fesErr errlist[] =
{
    { FES_SUCCESS,        "Success"},
    { FES_INPUT_ERROR,    "Missreading"},
    { FES_ACCESS_ERROR,   "Can't access grid"},
    { FES_GRIDS_CONSTANT, "The Definition of the grids is not constant"},
    { FES_NO_MEMORY,      "Out of memory"},
    { FES_NO_DATA,        "The tide is undefine"},
    { FES_MODE_ERROR,     "The requested mode to grids is not allowed"},
    { -1,                 NULL}
};

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
char* fesErr2String( const int err )
{
    int ix;

    for ( ix = 0; errlist[ix].code != -1; ix++ )
    {
	if ( err == errlist[ix].code )
	    return errlist[ix].reason;
    }

    return "Unknown error";
}
