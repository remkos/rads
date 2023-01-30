/* ######################################################################
 *  Main routines for the FES prediction software.
 * 
 *  File      : prediction.c
 *  Developer : CLS
 *  Version   : 1.1
 *  Date      : 6 oct 2004
 *  
 * ######################################################################
 */
#include "fes-int.h"

#define AAMU2	 0.069439968323
#define AANU2	-0.006104695053
#define AAL2	 0.077137765667
#define AAT2	 0.180480173707
#define AALDA2	 0.016503557465
#define BBMU2	 0.351535557706
#define BBNU2	 0.156878802427
#define BBL2	-0.051653455134
#define BBT2	-0.020101177502
#define BBLDA2	-0.013307812292
#define CCMU2	-0.046278307672
#define CCNU2	 0.006755704028
#define CCL2	 0.027869916824
#define CCT2	 0.008331518844
#define CCLDA2	 0.007753383202
#define AAP1	-0.238730166320
#define BBP1	 0.103860760343
#define CCP1	 0.289275485414





/*
// ///////////////////////////////////////////////////////////////////////////
// Computes square elevation.
//
// Parameters:
//	x:	value
//
// Return value:
//   value²
*/
static double sqr(const double x)
{
    return x * x;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Compute nodal corrections from SCHUREMAN (1958).
// indexes used in this routine are internal to the code
// and corresponds to the "original" ondes.dat file.
//
// Parameters:
//   formula
//   astro
*/
static void nodalA(const astronomicAngle* const astro,
		   wave* const w)
{
    int   ix;

    for ( ix = 0; ix < N_WAVES; ix++ )
    {
	switch ( w[ix].code )
	{
	case 1:
	case 27:
	case 65:
	case 66:
	case 67:
	case 69:
	    w[ix].f = sin(astro->iang) * sqr(cos(astro->iang / 2.0))/0.38;
	    break;
	case 2:
	case 12:
	case 13:
	case 60:
	case 71:
	case 72:
	    w[ix].f = 1.0;
	    break;
	case 3:
	    w[ix].f = sqrt(0.8965 * sqr(sin(2.0 * astro->iang)) +
		0.6001 * sin(2.0 * astro->iang) * cos(astro->nu) + 0.1006);
	    break;
	case 5:
	case 6:
	case 7:
	case 8:
	case 9:
	case 61:
	case 62:
	    w[ix].f = sqr(sqr(cos(astro->iang / 2.0))) / 0.9154;
	    break;
	case 11:
	    w[ix].f = sqr(sqr(cos(astro->iang / 2.0))) / 0.9154 * astro->x1ra;
	    break;
	case 14:
	    w[ix].f = sqrt(19.0444 * sqr(sqr(sin(astro->iang))) +
		2.7702 * sqr(sin(astro->iang)) * cos(2.0 * astro->nu) + 0.0981);
	    break;
	case 63:
	case 64:
	    w[ix].f = sin(astro->iang) * sin(astro->iang) / 0.1565;
	    break;
	case 68:
	case 70:
	case 73:
	case 74:
	    w[ix].f = sin(2.0 * astro->iang) / 0.7214;
	    break;
	case 75:
	    w[ix].f = sin(astro->iang) * sqr(sin(astro->iang / 2.0)) / 0.01640;
	    break;
	default:
	    assert(0);
	}
    }
}

/*
// ///////////////////////////////////////////////////////////////////////////
// Compute nodal corrections from SCHUREMAN (1958).
// indexes used in this routine are internal to the code
// and corresponds to the "original" ondes.dat file.
//
// Parameters:
//  w
//  astro
*/
static void nodalG(const astronomicAngle* const astro,
		   wave* const w)
{
    int                 ix;
    register double     u;
    register double     v0;

    for ( ix = 0; ix < N_WAVES; ix++ )
    {
	switch ( w[ix].code )
	{
	    /* O1 */
	case 1:
	    /*v0 = tt - 2.0 * s + hp + 90.0;
	    u = 2 * xi - nu;*/
	    w[ix].v0u = astro->tt - 2.0 * astro->s + astro->hp + 90.0 +
		2.0 * astro->xi - astro->nu;
	    break;
	    /* P1 */
	case 2:
	    /*v0 = tt - hp + 90.0;
	    u = 0.0;*/
	    w[ix].v0u=astro->tt - astro->hp + 90.0;
	    break;
	    /* K1 */
	case 3:
	    /*v0 = tt + hp - 90.0;
	    u = -nuprim;*/
	    w[ix].v0u = astro->tt + astro->hp - 90.0 - astro->nuprim;
	    break;
	    /* 2N2 */
	case 5:
	    /*v0 = 2.0 * tt - 4.0 * s + 2.0 * hp + 2.0 * p;
	    u = 2.0 * xi - 2.0 * nu;*/
	    w[ix].v0u =2.0 * astro->tt - 4.0 * astro->s + 2.0 * astro->hp +
		2.0 * astro->p + 2.0 * astro->xi - 2.0 * astro->nu;
	    break;
	    /* MU2 */
	case 6:
	    /*v0 = 2.0 * tt - 4.0 * s + 4.0 * hp;
	    u = 2.0 * xi - 2.0 * nu;*/
	    w[ix].v0u = 2.0 * astro->tt - 4.0 * astro->s + 4.0 * astro->hp +
		2.0 * astro->xi - 2.0 * astro->nu;
	    break;
	    /* N2 */
	case 7:
	    /*v0 = 2.0 * tt - 3.0 * s + 2.0 * hp + p;
	    u = 2.0 * xi - 2.0 * nu;*/
	    w[ix].v0u = 2.0 * astro->tt - 3.0 * astro->s + 2.0 * astro->hp +
		astro->p + 2.0 * astro->xi - 2.0 * astro->nu;
	    break;
	    /* Nu2 */
	case 8:
	    /*v0 = 2.0 * tt - 3.0 * s + 4.0 * hp - p;
	    u = 2.0 * xi - 2.0 * nu;*/
	    w[ix].v0u = 2.0 * astro->tt - 3.0 * astro->s + 4.0 * astro->hp -
		astro->p + 2.0 * astro->xi - 2.0 * astro->nu;
	    break;
	    /* M2 */
	case 9:
	    /*v0 = 2.0 * tt - 2.0 * s + 2.0 * hp;
	    u = 2.0 * xi - 2.0 * nu;*/
	    w[ix].v0u = 2.0 * astro->tt + - 2.0 * astro->s + 2.0 * astro->hp +
		2.0 * astro->xi - 2.0 * astro->nu;
	    break;
	    /* L2 */
	case 11:
	    /*v0 = 2.0 * tt - s + 2.0 * hp - p + 180.0;
	    u = 2.0 * xi - 2.0 * nu - r;*/
	    w[ix].v0u = 2.0 * astro->tt - astro->s + 2.0 * astro->hp -
		astro->p + 180.0 + 2.0 * astro->xi - 2.0 * astro->nu - astro->r;
	    break;
	    /* T2 */
	case 12:
	    /*v0 = 2.0 * tt - hp + p1;
	    u = 0.0;*/
	    w[ix].v0u = 2.0 * astro->tt - astro->hp + astro->p1;
	    break;
	    /* S2 */
	case 13:
	    /*v0 = 2.0 * tt;
	    u = 0.0;*/
	    w[ix].v0u = 2.0 * astro->tt;
	    break;
	    /* K2 */
	case 14:
	    /*v0 = 2.0 * tt + 2.0 * hp;
	    u = -2.0 * nusec;*/
	    w[ix].v0u = 2.0 * astro->tt + 2.0 * astro->hp - 2.0 * astro->nusec;
	    break;
	    /* Q1 */
	case 27:
	    /*v0 = tt - 3.0 * s + hp + p + 90.0;
	    u = 2.0 * xi - nu;*/
	    w[ix].v0u = astro->tt - 3.0 * astro->s + astro->hp + astro->p + 90.0 + 
		2.0 * astro->xi - astro->nu;
	    break;
	    /* Eps2 */
	case 60:
	    /*v0 = 2.0 * tt - 5 * s + 4 * hp + p;
	    u = 0.0;*/
	    w[ix].v0u = 2.0 * astro->tt - 5.0 * astro->s + 4.0 * astro->hp + 
		astro->p;
	    break;
	    /* Lambda2 */
	case 61:
	    v0 = 2.0 * astro->tt - astro->s + astro->p + 180.0;
	    u = 2.0 * astro->xi - 2.0 * astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* l21 */
	case 62:
	    v0 = 2.0 * astro->tt - astro->s + 2.0 * astro->hp - astro->p + 180.0;
	    u = -2.0 * astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* l22 */
	case 63:
	    v0 = 2.0 * astro->tt - astro->s + 2.0 * astro->hp + astro->p;
	    u = 2.0 * astro->xi - 2.0 * astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* Eta2 */
	case 64:
	    v0 = 2.0 * astro->tt + astro->s + 2.0 * astro->hp - astro->p;
	    u = -2.0 * astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* 2Q1 */
	case 65:
	    v0 = astro->tt - 4.0 * astro->s + astro->hp + 2.0 * astro->p + 90.0;
	    u = 2.0 * astro->xi - astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* Sigma1 */
	case 66:
	    v0 = astro->tt - 4.0 * astro->s + 3.0 * astro->hp + 90.0;
	    u = 2.0 * astro->xi - astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* Ro1 */
	case 67:
	    v0 = astro->tt - 3.0 * astro->s + 3.0 * astro->hp - astro->p + 90.0;
	    u = 2.0 * astro->xi - astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* M11 */
	case 68:
	    v0 = astro->tt - astro->s + astro->hp + astro->p - 90.0;
	    u = -astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* M12 */
	case 69:
	    v0 = astro->tt - astro->s + astro->hp - astro->p - 90.0;
	    u = 2.0 * astro->xi - astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* Ki1 */
	case 70:
	    v0 = astro->tt - astro->s + 3.0 * astro->hp - astro->p - 90.0;
	    u = -astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* Pi1 */
	case 71:
	    /*v0 = tt - 2.0 * hp + p1 + 90.0;
	    u = 0.0;*/
	    w[ix].v0u = astro->tt - 2.0 * astro->hp + astro->p1 + 90.0;
	    break;
	    /* Phi1 */
	case 72:
	    /*v0 = tt + 3.0 * hp - 90.0;
	    u = 0.0;*/
	    w[ix].v0u = astro->tt + 3.0 * astro->hp - 90.0;
	    break;
	    /* Teta1 */
	case 73:
	    v0 = astro->tt + astro->s - astro->hp + astro->p - 90.0;
	    u = -astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* J1 */
	case 74:
	    v0 = astro->tt + astro->s + astro->hp - astro->p - 90.0;
	    u = -astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	    /* OO1 */
	case 75:
	    v0 = astro->tt + 2.0 * astro->s + astro->hp - 90.0;
	    u = -2.0 * astro->xi - astro->nu;
	    w[ix].v0u = v0 + u;
	    break;
	default:
	    assert(0);
	}
	w[ix].v0u = fmod(w[ix].v0u, 360.00) * RAD;
    }
}


/*
// ///////////////////////////////////////////////////////////////////////////
// This program initialize some astronomic data useful for
// nodal corrections.
//   itj	Desired UTC time, in (decimal) Modified.
//   tt		Mean solar angle relative to Greenwich
//   hp		Mean solar longitude
//   s		Mean lunar longitude
//   p1		Longitude of solar perigee
//   p		Longitude of lunar perigee</param>
//   iang
//   xi
//   nu
//   x1ra
//   r
//   nuprim
//   nusec
*/
static void astronomics(const double  itj, astronomicAngle* angle)
{
    static const double ct0          = 180.0;
    static const double ct1          = 360.0 * 3.6525E+04;
    static const double cn0          = 259.1560563;
    static const double cn1          = -1934.1423972;
    static const double cs0          = 277.0256206;
    static const double cs1          = 481267.892;
    static const double ch0          = 280.1895015;
    static const double ch1          = 36000.76892;
    static const double cps0         = 281.2208568;
    static const double cps1         = 1.719175;
    static const double cp0          = 334.3837214;
    static const double cp1          = 4069.0322056;
    double              tgn2;
    double              at1;
    double              at2;
    double              u;
    double              tgi2;
    double              n;
    double              pp;

    /* tt mean solar angle relative to Greenwich */
    angle->tt = fmod(ct0 + ct1 * itj, 360.0) * RAD;

    /* hp longitude of ascending lunar node */
    n = fmod(cn0 + cn1 * itj, 360.0) * RAD;

    /* hp mean solar longitude */
    angle->hp = fmod(ch0 + ch1 * itj, 360.0) * RAD;

    /* s mean lunar longitude */
    angle->s = fmod(cs0 + cs1 * itj, 360.0) * RAD;

    /* p1 longitude of solar perigee */
    angle->p1 = fmod(cps0 + cps1 * itj, 360.0) * RAD;

    /* p longitude of lunar perigee */
    angle->p = fmod(cp0 + cp1 * itj, 360.0) * RAD;

    u = 9.13694997e-01 - 3.5692561e-02 * cos(n);

    angle->iang = acos(u);

    tgn2 = tan(n / 2.0);

    at1 = atan(1.01883 * tgn2);
    at2 = atan(6.4412e-01 * tgn2);

    angle->xi = -at1 - at2 + n;

    if ( n > M_PI )
    {
	angle->xi -= 2.0 * M_PI;
    }

    angle->nu = at1 - at2;

    /* for constituents l2,k1,k2 */
    tgi2 = tan(angle->iang / 2.0);

    pp = angle->p - angle->xi;

    angle->x1ra = sqrt(1.0 - 12.0 * sqr(tgi2) * cos(2.0 * pp) + 36.0 *
	pow(tgi2, 4));

    angle->r = atan(sin(2.0 * pp) / (1.0/(6.0 * sqr(tgi2)) - cos(2.0 * pp)));

    angle->nuprim = atan(sin(2.0 * (angle->iang)) * sin(angle->nu) /
	(sin(2.0 * (angle->iang)) * cos(angle->nu) + 3.347E-01));

    angle->nusec = 0.5 * atan(((sqr(sin(angle->iang))) * sin(2.0 * (angle->nu))) /
	(sqr(sin(angle->iang)) * cos(2.0 * (angle->nu))+ 7.27E-02));
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Compute nodal corrections.
//
// Parameters:
//  data:	Internal data
//  t0:		Desired UTC time, in (decimal) Modified.
*/
void initCorrections(const double t0,
		     wave* const w)
{
    astronomicAngle astro;

    astronomics(t0, &astro);

    nodalA(&astro, w);

    astro.tt            *= DEG;
    astro.hp            *= DEG;
    astro.s             *= DEG;
    astro.p1            *= DEG;
    astro.p             *= DEG;
    astro.xi            *= DEG;
    astro.nu            *= DEG;
    astro.r             *= DEG;
    astro.nuprim        *= DEG;
    astro.nusec         *= DEG;

    nodalG(&astro, w);

}



/*
// ///////////////////////////////////////////////////////////////////////////
// Init admittance coefficients
//
//   w:		  waves définition
//
*/
void initAdmittance(wave* const w)
{
    /* Q1		*/
    w[ 0].freq =       13.39866087990 * RAD;
    w[ 0].code =       27;
    /* O1		*/
    w[ 1].freq =       13.94303558000 * RAD;
    w[ 1].code =       1;
    /* K1		*/
    w[ 2].freq =       15.04106864000 * RAD;
    w[ 2].code =       3;
    /* 2N2	*/
    w[ 3].freq =       27.89535481990 * RAD;
    w[ 3].code =       5;
    /* N2		*/
    w[ 4].freq =       28.43972952010 * RAD;
    w[ 4].code =       7;
    /* M2		*/
    w[ 5].freq =       28.98410422000 * RAD;
    w[ 5].code =       9;
    /* K2		*/
    w[ 6].freq =       30.08213728000 * RAD;
    w[ 6].code =       14;
    /* S2		*/
    w[ 7].freq =       30.00000000000 * RAD;
    w[ 7].code =       13;
    /* P2		*/
    w[ 8].freq =       14.95893136000 * RAD;
    w[ 8].code =       2;
    /* Nu2		*/
    w[ 9].freq =       28.51258314000 * RAD;
    w[ 9].code =       8;
    /* Mu2		*/
    w[10].freq =       27.96820844000 * RAD;
    w[10].code =       6;
    /* L2		*/
    w[11].freq =       29.52847892000 * RAD;
    w[11].code =       11;
    /* T2		*/
    w[12].freq =       29.95893332010 * RAD;
    w[12].code =       12;
    /* Eps2	*/
    w[13].freq =       27.4238337 * RAD;
    w[13].code =       60;
    /* Lambda2	*/
    w[14].freq =       29.4556253 * RAD;
    w[14].code =       61;
    /* Eta2	*/
    w[15].freq =       30.6265120 * RAD;
    w[15].code =       64;
    /* 2Q1	*/
    w[16].freq =       12.8542862 * RAD;
    w[16].code =       65;
    /* Sigma1	*/
    w[17].freq =       12.9271398 * RAD;
    w[17].code =       66;
    /* Ro1	*/
    w[18].freq =       13.4715145 * RAD;
    w[18].code =       67;
    /* M11	*/
    w[19].freq =       14.4966939 * RAD;
    w[19].code =       68;
    /* M12	*/
    w[20].freq =       14.4874103 * RAD;
    w[20].code =       69;
    /* Ki1	*/
    w[21].freq =       14.5695476 * RAD;
    w[21].code =       70;
    /* Pi1	*/
    w[22].freq =       14.9178647 * RAD;
    w[22].code =       71;
    /* Phi1	*/
    w[23].freq =       15.1232059 * RAD;
    w[23].code =       72;
    /* teta1	*/
    w[24].freq =       15.5125897 * RAD;
    w[24].code =       73;
    /* J1		*/
    w[25].freq =       15.5854433 * RAD;
    w[25].code =       74;
    /* OO1	*/
    w[26].freq =       16.1391017 * RAD;
    w[26].code =       75;
}





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
{
    double              td;
    double              ph;
    double              shpn[4];
    double              zlp;

    /* Compute 4 principal mean longitudes in radians at time TD */
    td = ((ts + 33282.0) * 86400.0 - 4043174400.0) / 86400.0;

    shpn[0] = RAD * fmod(290.210 + (td * 13.17639650), 360.0);
    shpn[1] = RAD * fmod(280.120 + (td *  0.98564730), 360.0);
    shpn[2] = RAD * fmod(274.350 + (td *  0.11140410), 360.0);
    shpn[3] = RAD * fmod(343.510 + (td *  0.05295390), 360.0);

    /* Assemble long-period tide potential from 15 CTE terms > 1 mm.
    Nodal term is included but not the constant term. */
    zlp = 2.790 * cos(shpn[3]) -0.490 * cos(shpn[1] - (283.0 * RAD)) -3.100 *
	cos(2.00*shpn[1]);

    ph  = shpn[0];
    zlp = zlp - 0.670 * cos(ph - 2.0*shpn[1] + shpn[2]) -
	(3.520 - 0.460 * cos(shpn[3])) * cos(ph - shpn[2]);
    ph  = ph + shpn[0];

    zlp = zlp - 6.660 * cos(ph) - 2.760 * cos(ph + shpn[3]) - 0.260 *
	cos(ph + 2.0 * shpn[3]) - 0.580 * cos(ph - 2.0 * shpn[1]) - 0.290 *
	cos(ph - 2.0 * shpn[2]);
    ph  = ph + shpn[0];
    zlp = zlp - 1.270 * cos(ph - shpn[2]) - 0.530 *
	cos(ph - shpn[2] + shpn[3]) - 0.240 * cos(ph - 2.0 * shpn[1] + shpn[2]);

    /* Multiply by gamma_2 * sqrt(5/4 pi) * P20(lat) */
    *tlp = 0.4370 * zlp * ((1.50 * sqr(sin(lat * RAD))) - 0.50);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Compute the elapsed time since 01/01/1900 0:0:0.0 in julian centuries.
//
// Parameters:
//   date:	Modified Julian day (days since 1950-01-01 00:00:00.000 UTC)
//
// Return value:
//   The julian centuries.
*/
double julianCenturies(const double date)
{
    /* 18262.0 = number of day elapsed between 1950-01-01 00:00:00.000 UTC and
    1900-01-01 00:00:00.000 UTC */
    return(18262.0 + date) / 36525.0;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Compute admittance
//
// Parameters:
//   waves
*/
void admittance(const int computeP1,
		wave* const waves)
{
    register dComplex*  x;
    register dComplex*  y;
    register dComplex*  z;

    /* infer additional constituents by admittance DIURNALS
    (from Richard Ray perth2 program) */

    /* from Q1 and O1 (0-1) */

    x = &waves[0].cplx;
    y = &waves[1].cplx;
    z = &waves[2].cplx;

    /* 2Q1 */
    waves[16].cplx.re = 0.263 * x->re - 0.0252 * y->re;
    waves[16].cplx.im = 0.263 * x->im - 0.0252 * y->im;
    /* sigma1 */
    waves[17].cplx.re = 0.297 * x->re - 0.0264 * y->re;
    waves[17].cplx.im = 0.297 * x->im - 0.0264 * y->im;
    /* rho1 */
    waves[18].cplx.re = 0.164 * x->re + 0.0048 * y->re;
    waves[18].cplx.im = 0.164 * x->im + 0.0048 * y->im;

    /* from O1 and K1  (1-2) */

    /* M11 */
    waves[19].cplx.re =  0.0389 * y->re + 0.0282 * z->re;
    waves[19].cplx.im =  0.0389 * y->im + 0.0282 * z->im;
    /* M12 */
    waves[20].cplx.re =  0.0140 * y->re + 0.0101 * z->re;
    waves[20].cplx.im =  0.0140 * y->im + 0.0101 * z->im;
    /* chi1 */
    waves[21].cplx.re =  0.0064 * y->re + 0.0060 * z->re;
    waves[21].cplx.im =  0.0064 * y->im + 0.0060 * z->im;
    /* pi1 */
    waves[22].cplx.re =  0.0030 * y->re + 0.0171 * z->re;
    waves[22].cplx.im =  0.0030 * y->im + 0.0171 * z->im;
    /* phi1 */
    waves[23].cplx.re = -0.0015 * y->re + 0.0152 * z->re;
    waves[23].cplx.im = -0.0015 * y->im + 0.0152 * z->im;
    /* theta1 */
    waves[24].cplx.re = -0.0065 * y->re + 0.0155 * z->re;
    waves[24].cplx.im = -0.0065 * y->im + 0.0155 * z->im;
    /* J1 */
    waves[25].cplx.re = -0.0389 * y->re + 0.0836 * z->re;
    waves[25].cplx.im = -0.0389 * y->im + 0.0836 * z->im;
    /* OO1 */
    waves[26].cplx.re = -0.0431 * y->re + 0.0613 * z->re;
    waves[26].cplx.im = -0.0431 * y->im + 0.0613 * z->im;

    /* P1  from Grenoble admittance code */
    if ( computeP1 )
    {
	waves[8].cplx.re = AAP1 * x->re + BBP1 * y->re + CCP1 * z->re;
	waves[8].cplx.im = AAP1 * x->im + BBP1 * y->im + CCP1 * z->im;
    }

    /* SEMI-DIURNAL (from Grenoble to take advantage of 2N2) */

    /* from 2N2 -N2 (3-4) */

    x = &waves[3].cplx;
    y = &waves[4].cplx;

    /* eps2 */
    waves[13].cplx.re = 0.53285 * x->re - 0.03304 * y->re;
    waves[13].cplx.im = 0.53285 * x->im - 0.03304 * y->im;

    /* from M2 - K2 [5-6] */

    x = &waves[4].cplx;
    y = &waves[5].cplx;
    z = &waves[6].cplx;

    /* eta2 */
    waves[15].cplx.re = -0.0034925 * y->re + 0.0831707 * z->re;
    waves[15].cplx.im = -0.0034925 * y->im + 0.0831707 * z->im;

    /* from N2 -M2- K2 by spline admittances [see GRL 18[5]:845-848,1991] */

    /* mu2 */
    waves[10].cplx.re = AAMU2  * z->re + BBMU2  * x->re + CCMU2  * y->re;
    waves[10].cplx.im = AAMU2  * z->im + BBMU2  * x->im + CCMU2  * y->im;
    /* nu2 */
    waves[ 9].cplx.re = AANU2  * z->re + BBNU2  * x->re + CCNU2  * y->re;
    waves[ 9].cplx.im = AANU2  * z->im + BBNU2  * x->im + CCNU2  * y->im;
    /* lda2 */
    waves[14].cplx.re = AALDA2 * z->re + BBLDA2 * x->re + CCLDA2 * y->re;
    waves[14].cplx.im = AALDA2 * z->im + BBLDA2 * x->im + CCLDA2 * y->im;
    /* L2 */
    waves[11].cplx.re = AAL2   * z->re + BBL2   * x->re + CCL2   * y->re;
    waves[11].cplx.im = AAL2   * z->im + BBL2   * x->im + CCL2   * y->im;
    /* T2 */
    waves[12].cplx.re = AAT2   * z->re + BBT2   * x->re + CCT2   * y->re;
    waves[12].cplx.im = AAT2   * z->im + BBT2   * x->im + CCT2   * y->im;
}
