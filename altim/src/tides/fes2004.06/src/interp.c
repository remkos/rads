/*
 *  Bilinear interpolation.
 * 
 *  File      : interp.c
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.1
 *  Date      : 6 October 2004
 */

#include "fes-int.h"

/*
// ///////////////////////////////////////////////////////////////////////////
// Linear Weighting
//
// Parameters:
//	x:	abscissa where interpolation will occur
//	x1: 	abscissa corresponding to the first value
//	x2:	abscissa corresponding to the second value
//	w1	linear weight w1
//	w2:	linear weight w2
//
*/
static void linearWeighting(const double x,
			    const double x_1,
			    const double x_2,
			    double* w_1,
			    double* w_2)
{
    assert(x >= x_1 - EPSILON);
    assert(x <= x_2 + EPSILON);

    if ( EQUALS(x_1, x_2) || EQUALS(x, x_1) )
    {
	*w_1 = 1.0;
	*w_2 = 0.0;
    }
    else if ( EQUALS(x, x_2) )
    {
	*w_1 = 0.0;
	*w_2 = 1.0;
    }
    else
    {
	*w_1 = (x_2 - x) / (x_2 - x_1);
	*w_2 = (x - x_1) / (x_2 - x_1);
    }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Sum Weighting
//
// Parameters:
//	x_12
//	w_1
//	w_2
//	s
//	w
//
*/
static void sumWeighting(const double x12,
			 const double w12,
			 double* s,
			 double* w)
{
    if ( !EQUALS(x12, DV) )
    {
        *s += w12 * x12;
	*w += w12;
    }
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
*/
void bilinearInterp(const double x_1,
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
    double  w       = 0.0;
    double  s       = 0.0;
    double  w_X1;
    double  w_X2;
    double  w_Y1;
    double  w_Y2;

    linearWeighting(x, x_1, x_2, &w_X1, &w_X2);
    linearWeighting(y, y_1, y_2, &w_Y1, &w_Y2);
    
    sumWeighting(value_11, w_X1 * w_Y1, &s, &w);
    sumWeighting(value_12, w_X1 * w_Y2, &s, &w);
    sumWeighting(value_21, w_X2 * w_Y1, &s, &w);
    sumWeighting(value_22, w_X2 * w_Y2, &s, &w);
    
    *result = w == 0.0? DV: s / w;
}
