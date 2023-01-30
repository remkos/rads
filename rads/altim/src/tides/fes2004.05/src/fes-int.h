/*
 *  Include file for the FES prediction software.
 * 
 *  File      : fes-int.h
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.4
 *  Date      : 11 May 2005
 *  
 *  1.2: add M4, S1 (Ray), Mf, Mm, Mtm, MSqm
 *  1.3: remove conditional compilation
 *  1.4: correct computation of dynamic LP
 */
#ifndef _FESINT_H
#define _FESINT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#include <io.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <fcntl.h>
#include <math.h>
#include <assert.h>

#include "fes.h"

/* constants */
#define EPSILON			1.0E-9
#define DV			18446744073709551616.0
#define N_WAVES			33
#define DEPTH_TIDE		15
#define DEPTH_RADIAL		9
#define DEPTH			DEPTH_TIDE

#ifndef MAX_PATH
#define MAX_PATH		1024
#endif

#ifndef M_PI
#define M_PI			3.14159265358979323846
#endif

#define RAD			0.01745329251994329576
#define DEG			57.29577951308232087679

/* possible type of tide */
#define SP_TIDE			0
#define LP_TIDE			1


/* macros */
#define IN_GAP(a, b, c)		(((b) >= (a)) && ((b) <= (c)))
#define EQUALS(a, b)		(fabs(a - b) < EPSILON)


/*
 * structure representing a complex.
 */
typedef struct _fComplex
{
    float re;
    float im;
} fComplex;






/*
 * structure representing a complex.
 */
typedef struct _dComplex
{
    double  re;
    double  im;
} dComplex;





/*
 * structure representing the astronomical angles.
 */
typedef struct _astronomicAngle
{
    double  tt;
    double  hp;
    double  s;
    double  p1;
    double  p;
    double  iang;
    double  xi;
    double  nu;
    double  x1ra;
    double  r;
    double  nuprim;
    double  nusec;
} astronomicAngle;





/*
 * structure representing the waves definition
 */
typedef struct _wave
{
    int         code;
    int         type;
    double      freq;
    double      v0u;
    double      f;
    dComplex    cplx;
} wave;



/*
 * structure representing a grid.
 */
typedef struct _gridDsc
{
    int         littleEndian;
    int         tide;
    int         latSamples;
    int         lonSamples;
    int         depth;
    double      latMin;
    double      latMax;
    double      lonMin;
    double      lonMax;
    double      latStep;
    double      lonStep;
    double      undef;
    dComplex    wave[DEPTH];
    FILE**      handle;
    fComplex**  buffer;
} gridDsc;





/*
 *  structure for representing the FES prediction
 */
struct _fesData
{
    int         inGrid;
    int         isData;
    int         verbose;
    double      tNodal;
    double      southLat;
    double      northLat;
    double      westLon;
    double      eastLon;
    double      h[N_WAVES];
    dComplex    sw[DEPTH];
    dComplex    se[DEPTH];
    dComplex    nw[DEPTH];
    dComplex    ne[DEPTH];
    gridDsc     grid;
    wave        waves[N_WAVES];
};





/*
 * in interp.c
 */
void bilinearInterp( const double x_1, const double x_2, const double y_1,
		     const double y_2, const double value_11, const double value_21,
		     const double value_12, const double value_22, const double x,
		     const double y, double* result );

/*
 * in grid.c
 */
double getValue(const int idx, const double min, const double step);
int readHeader(FILE* handle, const int mode, const int iDepth, gridDsc* grid);
int interp(fesData* fes, const double lat, const double lon);

/*
 * in prediction.c
 */
void initAdmittance(wave* const w);
void initCorrections(const double t0, wave* const w);
double julianCenturies(const double date);
void admittance(wave* const waves);
void lpeqmt2(const double ts, const double lat, double* tlp);
#endif /* _FESINT_H */
