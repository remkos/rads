/* ######################################################################
 *  Main routines for the FES prediction software.
 * 
 *  Developer : CLS
 *  Version   : 1.1
 *  Date      : 6 oct 2004
 *  
 * ######################################################################
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
