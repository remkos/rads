/*
 *  Prediction routines.
 * 
 *  File      : prediction.c
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.3
 *  Date      : 11 May 2005
 *  
 *  1.2 : add M4, S1 (Ray), Mf, Mm, Mtm, MSqm
 *  1.3 : correct computation of LP
 */
#include "fes-int.h"

#define Q1	0
#define O1	1
#define K1	2
#define DN2	3
#define N2	4
#define M2	5
#define K2	6
#define S2	7
#define P1	8
#define M4	9
#define S1	10
#define MF	11
#define MM	12
#define MTM	13
#define MSQM	14
#define NU2	15
#define MU2	16
#define L2	17
#define T2	18
#define EPS2	19
#define LAMBDA2	20
#define ETA2	21
#define DQ1	22
#define SIGMA1	23
#define RO1	24
#define M11	25
#define M12	26
#define KI1	27
#define PI1	28
#define PHI1	29
#define TETA1	30
#define J1	31
#define OO1	32

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
	case 26:
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
	case 18:
	    w[ix].f = sqr(sqr(sqr(cos(astro->iang / 2.0))) / 0.9154);
	    break;
	case 40:
	case 42:
	case 45:
	    w[ix].f = sqr(sin(astro->iang)) / 0.1578;
	    break;
	case 38:
	    w[ix].f = (2.0 / 3.0 - sqr(sin(astro->iang))) / 0.5021;
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
	    w[ix].v0u = 2.0 * astro->tt - 2.0 * astro->s + 2.0 * astro->hp +
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
        /* M4 */
	case 18:
	    /*v0 = 4.0 * tt - 4.0 * s + 4.0 * hp;
	    u = 4.0 * xi - 4.0 * nu;*/
	    w[ix].v0u = 4.0 * astro->tt - 4.0 * astro->s + 4.0 * astro->hp +
		4.0 * astro->xi - 4.0 * astro->nu;
	    break;
	/* S1 */
	case 26:
	    /*v0 = 1.0 * tt;
	    u = 0.0;*/
	    w[ix].v0u = 1.0 * astro->tt;
	    break;
        /* Q1 */
	case 27:
	    /*v0 = tt - 3.0 * s + hp + p + 90.0;
	    u = 2.0 * xi - nu;*/
	    w[ix].v0u = astro->tt - 3.0 * astro->s + astro->hp + astro->p + 90.0 + 
		2.0 * astro->xi - astro->nu;
	    break;
	/* Mm */
	case 38:
	    /*v0 = s - p;
	    u = 0.0;*/
	    w[ix].v0u = astro->s - astro->p;
	    break;
	/* Mf */
	case 40:
	    /*v0 = 2.0 * s;
	    u = -2.0 * xi;*/
	    w[ix].v0u = 2.0 * astro->s - 2.0 * astro->xi;
	    break;
	/* Mtm */
	case 42:
	    /*v0 = 3.0 * s - p;
	    u = -2.0 * xi;*/
	    w[ix].v0u = 3.0 * astro->s - astro->p - 2.0 * astro->xi;
	    break;
	/* MSqm */
	case 45:
	    /*v0 = 4.0 * s - 2.0 * hp;
	    u = -2.0 * xi;*/
	    w[ix].v0u = 4.0 * astro->s - 2.0 * astro->hp - 2.0 * astro->xi;
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
    /* Q1 */
    w[Q1     ].freq =	13.39866087990 * RAD;
    w[Q1     ].code =	27;
    w[Q1     ].type =	SP_TIDE;
    /* O1 */
    w[O1     ].freq =	13.94303558000 * RAD;
    w[O1     ].code =	1;
    w[O1     ].type =	SP_TIDE;
    /* K1 */
    w[K1     ].freq =	15.04106864000 * RAD;
    w[K1     ].code =	3;
    w[K1     ].type =	SP_TIDE;
    /* 2N2 */
    w[DN2    ].freq =	27.89535481990 * RAD;
    w[DN2    ].code =	5;
    w[DN2    ].type =	SP_TIDE;
    /* N2 */
    w[N2     ].freq =	28.43972952010 * RAD;
    w[N2     ].code =	7;
    w[N2     ].type =	SP_TIDE;
    /* M2 */
    w[M2     ].freq =	28.98410422000 * RAD;
    w[M2     ].code =	9;
    w[M2     ].type =	SP_TIDE;
    /* K2 */
    w[K2     ].freq =	30.08213728000 * RAD;
    w[K2     ].code =	14;
    w[K2     ].type =	SP_TIDE;
    /* S2 */
    w[S2     ].freq =	30.00000000000 * RAD;
    w[S2     ].code =	13;
    w[S2     ].type =	SP_TIDE;
    /* P1 */
    w[P1     ].freq =	14.95893136000 * RAD;
    w[P1     ].code =	2;
    w[P1     ].type =	SP_TIDE;
    /* M4 */
    w[M4     ].freq =	57.96820840000 * RAD;
    w[M4     ].code =	18;
    w[M4     ].type =	SP_TIDE;
    /* S1 */
    w[S1     ].freq =	15.00000000000 * RAD;
    w[S1     ].code =	26;
    w[S1     ].type =	SP_TIDE;
    /* Mf */
    w[MF     ].freq =	1.098033100000 * RAD;
    w[MF     ].code =	40;
    w[MF     ].type =	LP_TIDE;
    /* Mm */
    w[MM     ].freq =	0.544374700000 * RAD;
    w[MM     ].code =	38;
    w[MM     ].type =	LP_TIDE;
    /* Mtm */
    w[MTM    ].freq =	1.642407800000 * RAD;
    w[MTM    ].code =	42;
    w[MTM    ].type =	LP_TIDE;
    /* MSqm */
    w[MSQM   ].freq =	2.113928800000 * RAD;
    w[MSQM   ].code =	45;
    w[MSQM   ].type =	LP_TIDE;
    /* Nu2 */
    w[NU2    ].freq =	28.51258314000 * RAD;
    w[NU2    ].code =	8;
    w[NU2    ].type =	SP_TIDE;
    /* Mu2 */
    w[MU2    ].freq =	27.96820844000 * RAD;
    w[MU2    ].code =	6;
    w[MU2    ].type =	SP_TIDE;
    /* L2 */
    w[L2     ].freq =	29.52847892000 * RAD;
    w[L2     ].code =	11;
    w[L2     ].type =	SP_TIDE;
    /* T2 */
    w[T2     ].freq =	29.95893332010 * RAD;
    w[T2     ].code =	12;
    w[T2     ].type =	SP_TIDE;
    /* Eps2 */
    w[EPS2   ].freq =	27.42383370000 * RAD;
    w[EPS2   ].code =	60;
    w[EPS2   ].type =	SP_TIDE;
    /* Lambda2 */
    w[LAMBDA2].freq =	29.45562530000 * RAD;
    w[LAMBDA2].code =	61;
    w[LAMBDA2].type =	SP_TIDE;
    /* Eta2 */
    w[ETA2   ].freq =	30.62651200000 * RAD;
    w[ETA2   ].code =	64;
    w[ETA2   ].type =	SP_TIDE;
    /* 2Q1 */
    w[DQ1    ].freq =	12.85428620000 * RAD;
    w[DQ1    ].code =	65;
    w[DQ1    ].type =	SP_TIDE;
    /* Sigma1 */
    w[SIGMA1 ].freq =	12.92713980000 * RAD;
    w[SIGMA1 ].code =	66;
    w[SIGMA1 ].type =	SP_TIDE;
    /* Ro1 */
    w[RO1    ].freq =	13.47151450000 * RAD;
    w[RO1    ].code =	67;
    w[RO1    ].type =	SP_TIDE;
    /* M11 */
    w[M11    ].freq =	14.49669390000 * RAD;
    w[M11    ].code =	68;
    w[M11    ].type =	SP_TIDE;
    /* M12 */
    w[M12    ].freq =	14.48741030000 * RAD;
    w[M12    ].code =	69;
    w[M12    ].type =	SP_TIDE;
    /* Ki1 */
    w[KI1    ].freq =	14.56954760000 * RAD;
    w[KI1    ].code =	70;
    w[KI1    ].type =	SP_TIDE;
    /* Pi1 */
    w[PI1    ].freq =	14.91786470000 * RAD;
    w[PI1].code =	71;
    w[PI1].type =	SP_TIDE;
    /* Phi1 */
    w[PHI1   ].freq =	15.12320590000 * RAD;
    w[PHI1   ].code =	72;
    w[PHI1   ].type =	SP_TIDE;
    /* teta1 */
    w[TETA1  ].freq =	15.51258970000 * RAD;
    w[TETA1  ].code =	73;
    w[TETA1  ].type =	SP_TIDE;
    /* J1 */
    w[J1     ].freq =	15.58544330000 * RAD;
    w[J1     ].code =	74;
    w[J1     ].type =	SP_TIDE;
    /* OO1 */
    w[OO1    ].freq =	16.13910170000 * RAD;
    w[OO1    ].code =	75;
    w[OO1    ].type =	SP_TIDE;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Computes the long-period equilibrium ocean tides.
// Processing logic -
//   Fifteen tidal spectral lines from the Cartwright-Tayler-Edden
//   tables are summed over to compute the long-period tide.
//   Waves Mm, Mf and Mtm are removed (CLS 2004).
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
void lpeqmt2(const double ts, const double lat, double* tlp)
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
    /* Mm Removed
    zlp = zlp - 0.670 * cos(ph - 2.0*shpn[1] + shpn[2]) -
	(3.520 - 0.460 * cos(shpn[3])) * cos(ph - shpn[2]); */
    zlp = zlp - 0.670 * cos(ph - 2.0*shpn[1] + shpn[2]) +
	(0.460 * cos(shpn[3])) * cos(ph - shpn[2]);
    ph  = ph + shpn[0];

    /* Mf removed
    zlp = zlp - 6.660 * cos(ph) - 2.760 * cos(ph + shpn[3]) - 0.260 *
	cos(ph + 2.0 * shpn[3]) - 0.580 * cos(ph - 2.0 * shpn[1]) - 0.290 *
	cos(ph - 2.0 * shpn[2]); */
    zlp = zlp - 2.760 * cos(ph + shpn[3]) - 0.260 *
	cos(ph + 2.0 * shpn[3]) - 0.580 * cos(ph - 2.0 * shpn[1]) - 0.290 *
	cos(ph - 2.0 * shpn[2]);
    ph  = ph + shpn[0];
    /* Mtm removed
    zlp = zlp - 1.270 * cos(ph - shpn[2]) - 0.530 *
	cos(ph - shpn[2] + shpn[3]) - 0.240 * cos(ph - 2.0 * shpn[1] + shpn[2]); */
    zlp = zlp - 0.530 *
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
void admittance(wave* const waves)
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
    waves[DQ1   ].cplx.re = 0.263 * x->re - 0.0252 * y->re;
    waves[DQ1   ].cplx.im = 0.263 * x->im - 0.0252 * y->im;
    /* sigma1 */
    waves[SIGMA1].cplx.re = 0.297 * x->re - 0.0264 * y->re;
    waves[SIGMA1].cplx.im = 0.297 * x->im - 0.0264 * y->im;
    /* ro1 */
    waves[RO1   ].cplx.re = 0.164 * x->re + 0.0048 * y->re;
    waves[RO1   ].cplx.im = 0.164 * x->im + 0.0048 * y->im;

    /* from O1 and K1  (1-2) */

    /* M11 */
    waves[M11  ].cplx.re =  0.0389 * y->re + 0.0282 * z->re;
    waves[M11  ].cplx.im =  0.0389 * y->im + 0.0282 * z->im;
    /* M12 */
    waves[M12  ].cplx.re =  0.0140 * y->re + 0.0101 * z->re;
    waves[M12  ].cplx.im =  0.0140 * y->im + 0.0101 * z->im;
    /* KI1 */
    waves[KI1  ].cplx.re =  0.0064 * y->re + 0.0060 * z->re;
    waves[KI1  ].cplx.im =  0.0064 * y->im + 0.0060 * z->im;
    /* pi1 */
    waves[PI1  ].cplx.re =  0.0030 * y->re + 0.0171 * z->re;
    waves[PI1  ].cplx.im =  0.0030 * y->im + 0.0171 * z->im;
    /* phi1 */
    waves[PHI1 ].cplx.re = -0.0015 * y->re + 0.0152 * z->re;
    waves[PHI1 ].cplx.im = -0.0015 * y->im + 0.0152 * z->im;
    /* teta1 */
    waves[TETA1].cplx.re = -0.0065 * y->re + 0.0155 * z->re;
    waves[TETA1].cplx.im = -0.0065 * y->im + 0.0155 * z->im;
    /* J1 */
    waves[J1   ].cplx.re = -0.0389 * y->re + 0.0836 * z->re;
    waves[J1   ].cplx.im = -0.0389 * y->im + 0.0836 * z->im;
    /* OO1 */
    waves[OO1  ].cplx.re = -0.0431 * y->re + 0.0613 * z->re;
    waves[OO1  ].cplx.im = -0.0431 * y->im + 0.0613 * z->im;

    /* SEMI-DIURNAL (from Grenoble to take advantage of 2N2) */

    /* from 2N2 -N2 (3-4) */

    x = &waves[DN2].cplx;
    y = &waves[N2].cplx;

    /* eps2 */
    waves[EPS2].cplx.re = 0.53285 * x->re - 0.03304 * y->re;
    waves[EPS2].cplx.im = 0.53285 * x->im - 0.03304 * y->im;

    /* from M2 - K2 [5-6] */

    x = &waves[N2].cplx;
    y = &waves[M2].cplx;
    z = &waves[K2].cplx;

    /* eta2 */
    waves[ETA2].cplx.re = -0.0034925 * y->re + 0.0831707 * z->re;
    waves[ETA2].cplx.im = -0.0034925 * y->im + 0.0831707 * z->im;

    /* from N2 -M2- K2 by spline admittances [see GRL 18[5]:845-848,1991] */

    /* mu2 */
    waves[MU2    ].cplx.re = AAMU2  * z->re + BBMU2  * x->re + CCMU2  * y->re;
    waves[MU2    ].cplx.im = AAMU2  * z->im + BBMU2  * x->im + CCMU2  * y->im;
    /* nu2 */
    waves[NU2    ].cplx.re = AANU2  * z->re + BBNU2  * x->re + CCNU2  * y->re;
    waves[NU2    ].cplx.im = AANU2  * z->im + BBNU2  * x->im + CCNU2  * y->im;
    /* lambda2 */
    waves[LAMBDA2].cplx.re = AALDA2 * z->re + BBLDA2 * x->re + CCLDA2 * y->re;
    waves[LAMBDA2].cplx.im = AALDA2 * z->im + BBLDA2 * x->im + CCLDA2 * y->im;
    /* L2 */
    waves[L2     ].cplx.re = AAL2   * z->re + BBL2   * x->re + CCL2   * y->re;
    waves[L2     ].cplx.im = AAL2   * z->im + BBL2   * x->im + CCL2   * y->im;
    /* T2 */
    waves[T2     ].cplx.re = AAT2   * z->re + BBT2   * x->re + CCT2   * y->re;
    waves[T2     ].cplx.im = AAT2   * z->im + BBT2   * x->im + CCT2   * y->im;
}
