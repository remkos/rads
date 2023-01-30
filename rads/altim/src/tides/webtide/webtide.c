/*WEBTIDE -- Compute tide from BIO/OSD tide model
*+
      INTEGER FUNCTION WEBTIDEINIT (DIRNAME, EXT, ID)
      CHARACTER*(*) DIRNAME, EXT
      INTEGER ID

      INTEGER FUNCTION WEBTIDE (ID, UTC, LAT, LON, TIDE, TIDE2)
      INTEGER ID
      REAL*8 UTC, LAT, LON, TIDE, TIDE2

      SUBROUTINE WEBTIDEFREE (ID)
      INTEGER ID

* The Bedford Insitute of Ocenography (BIO) Ocean Science Division (OSD) of the
* Department of Fisheries and Oceans Canada provides various local tide models
* for the coasts of the North American Continent. These models can be found at:
* http://www.mar.dfo-mpo.gc.ca/science/ocean/coastal_hydrodynamics/WebTide/webtide.html
*
* The WEBTIDE routine has to be initialised using WEBTIDEINIT. Use the name of
* the directory and the preferred file extension as arguments. The directory
* DIRNAME should contain the following files:
* - tidecor.cfg, which contains the names of the node, elements and IOS file
* - The files mentioned in tidecor.cfg
* - constituents.txt, which contains the names of all constituents ending with "NONE"
* - Files CONST.barotropic.EXT, where CONST is the name of any consituent and EXT is
*   the preferred extension.
* The extension EXT can be:
* - s2c : tidal elevation
* - v2c : tidal velocities
* - v2r : relative tide (add Z0 component)
* WEBTIDEINIT returns an ID number (between 1 and 9) that needs to be used in
* subsequent calls to WEBTIDE. This provides the possibility to initialize
* and evaluate up to 9 models at a time. When no more slots are available ID = 0
* is returned.
*
* WEBTIDE can be called to determine the tidal elevation or velocity at any time
* (UTC) and location (LAT, LON). If the tidal model called for is a velocity model
* (v2c files instead of s2c files), both TIDE and TIDE2 will be filled.
* When the location is outside the model (that is outside any of the finite
* elements), TIDE (and TIDE2) will be filled with NaN (Not-A-Number).
*
* After the last call to WEBTIDE use WEBTIDEFREE to free the memory occupied by
* the tide model data. Again, use the ID number to indicate which model to free,
* or use 0 to free all models simultaneously.
*
* Arguments:
*   DIRNAME: Path name of the directory containing the requested WEBTIDE tide
*            model
*   EXT    : Extension of the constituent file names (s2c, v2c, v2r)
*   ID     : Identifier for WebTide model returned by WEBTIDEINIT
*   UTC    : Tide in UTC seconds since 1.0 January 1985
*   LAT    : Latitude in degrees
*   LON    : Longitude in degrees
*   TIDE   : Tide in meters
* In case the tide model is complex (.v2c files are used):
*   TIDE   : X-component of tidal velocity (m/s)
*   TIDE2  : Y-component of tidal velocity (m/s)
*
* Return codes:
*   WEBTIDEINIT: 0 = Scalar tide, 1 = Complex tide
*   WEBTIDE    : Returns index number of the element in which the location
*                resides. Returns -1 if outside model region.
*
* This code can also be called from C programs using the functions WebTideInit,
* WebTide and WebTideFree.
*=
* The code was originally written by M. Foreman as an adaptation of "tidecor" code by
* Patrick Roussel. Versions 2.0 and later of that code were created by Jason Chaffey.
*
* $Log: webtide.c,v $
* Revision 1.3  2011/06/16 18:23:27  rads
* - Allow multiple instances of WebTide simultaneously using the W struct in the C
*   routines and the ID identifier in the Fortran routines.
* - Major code tidying.
*
* Revision 1.2  2007/03/05 20:19:02  rads
* - Make NaN using make_dnan for porting to Linux
*
* Revision 1.1  2007/03/05 15:17:09  rads
* Initial revision of webtide code
*
* The subroutines below have been adapted by Remko Scharroo (Altimetrics LLC) from the
* program tidecor2.42.c provided at the BIO/OSD web site. The changes to the original
* software include:
*
* - Created the routines WebTideInit, WebTide and WebTideFree from the original main
*   program in tidecor2.42.c.
* - Allow running multiple instances simultaneously through the W struct in the C
*   routines or the id identifier in the Fortran routines.
* - Created the Fortran wrappers: WEBTIDEINIT, WEBTIDE and WEBTIDEFREE.
* - Added much more usage documentation.
* - Removed the minimum tide (can be reinserted by defining MIN_TIDE).
* - Sped up of the code by searching only those elements that can possibly contain
*   the closest node by storing the minimum and maximum element that contains a node.
* - No more searching of all elements after testing the closest elements fails (unless
*   SEARCH_ALL is defined). This makes barely any difference. A test of 40000 coastal
*   points only had 15 points found using this piece of code.
* - More acceleration by checking on the model boundaries in advance of the node search.
* - Return NaN (Not-A-Number) when location is outside any finite element (instead of
*   returning the tide at the closest point).
* - Fixed memory leaks (the original code added 24 bytes in phi2d on every iteration).
* - Do not use a special configuration file, but the files available in the model
*   directory.
* - A lot more tidying up of the code.
*-----------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "webtide.h"

/* Define the value for 2 * pi */
#define TWO_PI	6.28318530717958647692
#define DEG2RAD	(TWO_PI / 360.0)

/* To make a double NaN */
#define make_dnan(x) (((unsigned int *) &x)[0] = 0xffffffff, ((unsigned int *) &x)[1] = 0xffffffff)

/* These define the indices for the 5 constituents used in this program
   (Perhaps someday users will be allowed to choose which and how many to use) */
#define m2  0
#define s2  1
#define n2  2
#define k1  3
#define o1  4

/* Global constants: */
const int name_len = 6;		/* Max length of constituent names */
const int maxcons = 45;		/* Initial max # of main constituents */
const int maxshll = 110;	/* Initial max # of shallow water constituents */

/* Stack of memory locations for Fortran interface */
#define maxmodels 9
struct WebTideInfo X[maxmodels];
static int first = 1;

/* Prototype declarations: */
struct satptr *add_sat (char *satin);
void astro_angles (double d1, double *h, double *pp, double *s, double *p,
	double *np, double *dh, double *dpp, double *ds, double *dp, double *dnp);
int basis2d (struct WebTideInfo *W, double ptx, double pty, double *basis, int *in);
char *caseup (char *lower);
int closestnode (struct WebTideInfo *W, double ptx, double pty);
double dis_sq (double lng1, double lat1, double lng2, double lat2, double coslat);
double dis_man (double lng1, double lat1, double lng2, double lat2, double coslat);
int openvuf (struct WebTideInfo *W, char *iosfile);
void phi2d (double *xloc, double *yloc, double ptx, double pty, double *basis);
int raybound (double *xd, double *yd, double ptx, double pty);
void readnodes (struct WebTideInfo *W, char *nfile, char *efile);
void setvuf (struct WebTideInfo *W, long int kd, double xlat);
double TideP (struct WebTideInfo *W, long int kd, double dthr, double latitude, double *ampl, double *phase);
int vuf (struct WebTideInfo *W, char *inname);

#include <sysdep.h>

#define stringlength 1024

#ifdef CAPITALS
#define webtide WEBTIDE
#define webtideinit WEBTIDEINIT
#define webtidefree WEBTIDEFREE
#endif

#ifdef UNDERSCOREAFTER
#define webtide webtide_
#define webtideinit webtideinit_
#define webtidefree webtidefree_
#endif

#ifdef UNDERSCOREBEFORE
#define webtide _webtide
#define webtideinit _webtideinit
#define webtidefree _webtidefree
#endif

/*****************************************************************************/
/* These are Fortran Wrappers */

long int webtideinit (char *dirname, char *ext, int *id, int len1, int len2)
{
	char temp1[stringlength], temp2[4];        /* Scratch string...could malloc () */
	int len, i;

	/* Truncate Fortran string and end with null-char */
	strncpy (temp1, dirname, len1);
	len = len1 - 1;
	while (temp1[len] == ' ' && len >= 0) len--;
	temp1[len+1] = '\0';
	len = len2 - 1;
	strncpy (temp2, ext, len2);
	while (temp2[len] == ' ' && len >= 0) len--;
	temp2[len+1] = '\0';

	/* On first call, init */
	if (first) {
		for (i = 0; i < maxmodels; i++) X[i].nconst = 0;
		first = 0;
	}

	/* Look for free memory slot */
	*id = 0;
	for (i = 0; i < maxmodels; i++) {
		if (X[i].nconst > 0) continue;
		*id = i+1;
		return (WebTideInit (temp1, temp2, &X[i]));
	}

	fprintf (stderr, "webtideinit: no more memory slots available\n");
	return (0);
}

long int webtide (int *id, double *time, double *latitude, double *longitude, double *tide, double *tide2)
{
	return (WebTide (&X[*id-1], time, latitude, longitude, tide, tide2));
}

void webtidefree (int *id)
{
	int i, i0, i1;

	if (*id > 0 && *id <= maxmodels)
		i0 = *id-1, i1 = *id;		/* Free only one model */
	else
		i0 = 0, i1 = maxmodels;		/* Free all models */

	for (i = i0; i < i1; i++) {
		WebTideFree (&X[i]);
		X[i].nconst = 0;
	}
}

/*****************************************************************************/
long int WebTideInit (char *dirname, char *ext, struct WebTideInfo *W)
{
	char *dummy, *nodefile, *elemfile, *iosfile, *ext1;
	double *amp, *amp2, *phase, *phase2;
#ifdef MIN_TIDE
	double *min_tide, *min_tide2;
#endif
	FILE *fd, **indatfile;
	int i, j, garbage, nread;

/* Check the data file format */

	ext1 = (char *)malloc (5 * sizeof (char));
	ext1 = caseup (ext);
	if (!strcmp (ext1, "V2C"))
		W->cplx = 1;
	else if (!strcmp (ext1, "S2C"))
		W->cplx = 0;
	else if (!strcmp (ext1, "V2R"))
		W->cplx = -1;
	else {
		fprintf (stderr, "Error. Unrecognized data file format (%s).\n", ext1);
		exit (-12);
	}

/* Make space for file names */

	dummy    = (char *)malloc (180 * sizeof (char));
	nodefile = (char *)malloc (180 * sizeof (char));
	elemfile = (char *)malloc (180 * sizeof (char));
	iosfile  = (char *)malloc (180 * sizeof (char));

/* Read in file names of nodes, elements and IOS table from configuration file */

	sprintf (dummy, "%s/tidecor.cfg", dirname);
	if ( (fd = fopen (dummy, "r")) == NULL) {
		fprintf (stderr,"Error in opening the configuration file (%s) for WebTide\n", dummy);
		fprintf (stderr, "Cannot continue!\n\n");
		exit (9);
	}

/* Now read other filenames */

	j = fscanf (fd, "%s", dummy);
	if (j == EOF || j <= 0) {
		fprintf (stderr,"Error in reading the node filename.\n");
		fprintf (stderr, "Cannot continue!\n\n");
		exit (9);
	}
	sprintf (nodefile, "%s/%s", dirname, dummy);
	j = fscanf (fd, "%s", dummy);
	if (j == EOF || j <= 0) {
		fprintf (stderr,"Error in reading the element filename.\n");
		fprintf (stderr, "Cannot continue!\n\n");
		exit (9);
	}
	sprintf (elemfile, "%s/%s", dirname, dummy);
	j = fscanf (fd, "%s", dummy);
	if (j == EOF || j <= 0) {
		fprintf (stderr,"Error in reading the IOS table filename.\n");
		fprintf (stderr, "Cannot continue!\n\n");
		exit (9);
	}
	sprintf (iosfile, "%s/%s", dirname, dummy);
	fclose (fd);

	/* Read in the mesh */

	readnodes (W, nodefile, elemfile);

	/* Open the constituents file and determine the number of constants */

	sprintf (dummy, "%s/constituents.txt", dirname);
	if ( (fd = fopen (dummy, "r")) == NULL) {
		fprintf (stderr,"Error in opening the consituents file (%s) for WebTide\n", dummy);
		fprintf (stderr, "Cannot continue!\n\n");
		exit (9);
	}
	W->nconst = -1;
	while (!feof (fd)) {
		W->nconst++;
		fscanf (fd, "%s", dummy);
		if (!strcmp (dummy, "NONE")) break;
	}
	rewind (fd);

	/* Read all the constituent data files */

	W->constnames = (char **)malloc (W->nconst * sizeof (char *));
	for (i = 0; i < W->nconst; i++) {
		W->constnames[i]=(char *)malloc (name_len * sizeof (char));
	}
	indatfile = (FILE **)malloc (W->nconst * sizeof (FILE *));

	for (i = 0; i < W->nconst; i++) {
		j = fscanf (fd, "%s", W->constnames[i]);
		if (j == EOF || j <= 0) {
			fprintf (stderr,"Error in reading constituent #%d name.\n", i);
			fprintf (stderr, "Cannot continue!\n\n");
			exit (9);
		}
		W->constnames[i] = caseup (W->constnames[i]);
		sprintf (dummy, "%s/%s.barotropic.%s", dirname, W->constnames[i], ext);
		if ( (indatfile[i] = fopen (dummy, "r")) == NULL) {
			fprintf (stderr, "Error in opening constituent datafile (%s)\n", dummy);
			fprintf (stderr, "Cannot continue!\n\n");
			exit (11);
		}
		fgets (dummy, 180, indatfile[i]);
		fgets (dummy, 180, indatfile[i]);
		if (strcmp (W->constnames[i], "Z0") != 0) fgets (dummy, 180, indatfile[i]);
	}
	fclose (fd);

/* Allocate memory */

	W->cons = (struct mainptr**)malloc (maxcons * sizeof (struct mainptr*));
	W->shall = (struct shallptr**)malloc (maxshll * sizeof (struct shallptr*));
	W->amp_base = (double *)malloc (W->numnodes * W->nconst * sizeof (double));
	W->phase_base = (double *)malloc (W->numnodes * W->nconst * sizeof (double));
#ifdef MIN_TIDE
	W->min_tide_base = (double *)malloc (W->numnodes * sizeof (double));
#endif
	if (W->cplx > 0) {
		W->amp2_base = (double *)malloc (W->numnodes * W->nconst * sizeof (double));
		W->phase2_base = (double *)malloc (W->numnodes * W->nconst * sizeof (double));
#ifdef MIN_TIDE
		W->min_tide2_base = (double *)malloc (W->numnodes * sizeof (double));
#endif
	}
	else {
		W->amp2_base = W->phase2_base = NULL;
#ifdef MIN_TIDE
		W->min_tide2_base = NULL;
#endif
	}

/* Load the model tidal data */

	amp = W->amp_base;
	amp2 = W->amp2_base;
	phase = W->phase_base;
	phase2 = W->phase2_base;
	for (j = 0; j < W->numnodes; j++) {
		for (i = 0; i < W->nconst; i++) {
			if (W->cplx == 0) {
				if ( (nread = fscanf (indatfile[i], "%d %lf %lf", &garbage, amp++, phase++)) != 3) {
					fprintf (stderr, "Error reading phase and amplitude file!!\n");
					exit (1);
				}
			}
			else if (strcmp (W->constnames[i], "Z0") == 0) {
				if ( (nread = fscanf (indatfile[i], "%d %lf %lf", &garbage, amp++, amp2++)) != 3) {
					fprintf (stderr, "Error reading mean current file %d!!\n", i);
					exit (1);
				}
				*phase = 0.0; phase++;
				*phase2 = 0.0; phase2++;
			}
			else {
				if ( (nread = fscanf (indatfile[i], "%d %lf %lf %lf %lf", &garbage, amp++, phase++, amp2++, phase2++)) != 5) {
					fprintf (stderr, "Error reading vel phase and amplitude file %d!!\n", i);
					exit (1);
				}
			}
		}
	}

#ifdef MIN_TIDE
/* Calculate the minimum tides for each node */

	amp = W->amp_base;
	amp2 = W->amp2_base;
	min_tide = W->min_tide_base;
	min_tide2 = W->min_tide2_base;
	for (j = 0; j < W->numnodes; j++) {
		*min_tide = 0.0;
		for (i = 0; i < W->nconst; i++) *min_tide += *amp++;
		min_tide++;
		if (W->cplx > 0) {
			*min_tide2 = 0.0;
			for (i = 0; i < W->nconst; i++) *min_tide2 += *amp2++;
			min_tide2++;
		}
	}
#endif

/* Read in the constituent data */
	openvuf (W, iosfile);

/* Close the model tidal data files */
	for (i = 0; i < W->nconst; i++) fclose (indatfile[i]);
	free (nodefile);
	free (elemfile);
	free (iosfile);
	free (dummy);
	free (ext1);
	free (indatfile);

	return (W->cplx);
}

/*****************************************************************************/

long WebTide (struct WebTideInfo *W, double *time, double *latitude, double *longitude, double *tide, double *tide2)
{
	double lon, dt, elem_res[3], elem_res2[3], basis[3], *amp, *phase;
	int kd, i, elem, *in;

	lon = *longitude;
	if (lon > 180.0) lon -= 360.0;
	elem = basis2d (W, lon, *latitude, basis, in);

	if (elem < 0) {
/* No element was found that contained this position.
   Return NAN (Not-A-Number) */
		make_dnan (*tide);
		if (W->cplx > 0) *tide2 = *tide;
		return (elem);
	}

/* Split UTC85 into days since 1900.0 (kd) and fraction of day in hours (dt) */
	dt = *time / 86400.0 + 31047.0;
	kd = (long int)floor (dt);
	dt = (dt - (double)kd) * 24.0;

/* If the point is inside an element, calculate the tidal correction for
   each node of the element and interpolate to get the tidal correction
   at the new position. */

	for (i = 0; i <= 2; i++) {
		in = W->in_base + i + (3 * elem);
		amp = W->amp_base + (*in) * W->nconst;
		phase = W->phase_base + (*in) * W->nconst;
		elem_res[i] = TideP (W, kd, dt, W->ydata[*in], amp, phase);
#ifdef MIN_TIDE
		min_tide = W->min_tide_base + *in;
		elem_res[i] += *min_tide;
#endif
		if (W->cplx > 0) {
			amp = W->amp2_base + (*in) * W->nconst;
			phase = W->phase2_base + (*in) * W->nconst;
			elem_res2[i] = TideP (W, kd, dt, W->ydata[*in], amp, phase);
#ifdef MIN_TIDE
			amp = W->min_tide2_base + *in;
			elem_res2[i] += *amp;
#endif
		}
	}
	*tide = elem_res[0] * basis[0] + elem_res[1] * basis[1] + elem_res[2] * basis[2];
	if (W->cplx > 0) *tide2 = elem_res2[0] * basis[0] + elem_res2[1] * basis[1] + elem_res2[2] * basis[2];
	return (elem);
}

/*****************************************************************************/
void WebTideFree (struct WebTideInfo *W)
{
	int i, j;
	free (W->xdata);
	free (W->ydata);
	free (W->minelem);
	free (W->maxelem);
	free (W->in_base);
	free (W->amp_base);
	free (W->phase_base);
#ifdef MIN_TIDE
	free (W->min_tide_base);
#endif
	if (W->cplx > 0) {
		free (W->amp2_base);
		free (W->phase2_base);
#ifdef MIN_TIDE
		free (W->min_tide2_base);
#endif
	}

	for (i = 0; i < W->nconst; i++) {
		free (W->constnames[i]);
	}
	free (W->constnames);

	for (i = 0; i < W->nmain; i++) {
		for (j = 0; j < W->cons[i]->nsats; j++) {
			if (W->cons[i]->sats[j]) free (W->cons[i]->sats[j]);
		}
		if (W->cons[i]->nsats > 0) free (W->cons[i]->sats);
		free (W->cons[i]->name);
		free (W->cons[i]);
	}
	free (W->cons);

	for (i = 0; i < W->nshall; i++) {
		for (j = 0; j < W->shall[i]->numcon ; j++) {
			free (W->shall[i]->shcon[j]->name);
			free (W->shall[i]->shcon[j]);
		}
		free (W->shall[i]->shcon);
		free (W->shall[i]->name);
		free (W->shall[i]);
	}
	free (W->shall);
}


/*****************************************************************************/
struct satptr *add_sat (char *satin)
/* Gets satellite info from input line, puts the info into a new
   satellite structure and returns the structure. */
{
	char *dummy;
	int d1, d2, d3;
	double phse, rat;
	struct satptr *newsat;

	sscanf (satin, "%d %d %d %lf %lf", &d1,&d2,&d3, &phse, &rat);

	newsat = (struct satptr *)malloc (sizeof (struct satptr));

	newsat->deld[0] = d1;
	newsat->deld[1] = d2;
	newsat->deld[2] = d3;
	newsat->phase = phse;
	newsat->ratio = rat;

	if ( (dummy = strchr (satin, 82)) == NULL)
		newsat->corr = 0;
	else if (strncmp (dummy, "R1", 2) == 0)
		newsat->corr = 1;
	else if (strncmp (dummy, "R2", 2) == 0)
		newsat->corr = 2;
	else
		newsat->corr = 0;

	return newsat;
}

/*****************************************************************************/
void astro_angles (double d1, double *h, double *pp, double *s, double *p,
		double *np, double *dh, double *dpp, double *ds, double *dp, double *dnp)
/* Calculates the following ephermides of the sun and moon:
   h  = mean longitude of the sun;
   pp = mean longitude of the solar perigee;
   s  = mean longitude of the moon;
   p  = mean longitude of the lunar perigee; and
   np = negative of the longitude of the mean ascending node.
   Also calculates their rate of change.
   Units are cycles (cycles / 365 days for rates of change).
   These formulae were taken from pp 98 and 107 of the Explanatory
   Supplement to the Astronomical Ephermeris and the American
   Ephermis and Nautical Almanac (1961)
*/
{
	double dum, d12, d2, d22, d23, f, f2;

	d12 = d1 * d1;
	d2 = d1 * 1.0E-04;
	d22 = d2 * d2;
	d23 = d22 * d2;
	f = 360.0;
	f2 = f / 365.0;

	*h = (2.79696678E+02 + d1 * 9.856473354E-01 + d22 * 2.267E-05) / f;
	*h = modf (*h, &dum);

	*pp = (2.81220833E+02 + d1 * 4.70684E-05 + d22 * 3.39E-05 + d23 * 7.0E-08) / f;
	*pp = modf (*pp, &dum);

	*s = (2.70434164E+02 + d1 * 1.31763965268E+01 - d22 * 8.5E-05 + d23 * 3.9E-08) / f;
	*s = modf (*s, &dum);

	*p = (3.34329556E+02 + d1 * 1.114040803E-01 - d22 * 7.739E-04 - d23 * 2.6E-07) / f;
	*p = modf (*p, &dum);

	*np = (-2.59183275E+02 + d1 * 5.29539222E-02 - d22 * 1.557E-04 - d23 * 5.0E-08) / f;
	*np = modf (*np, &dum);

	*dh = (9.856473354E-01 + d1 * 2.267E-05 * 2.0E-08) / f2;

	*dpp = (4.70684E-05 + d1 * 3.39E-05 * 2.0E-08 + d12 * 7.0E-08 * 3.0E-12) / f2;

	*ds = (1.31763965268E+01 - d1 * 8.5E-05 * 2.0E-08 + d12 * 3.9E-08 *3.0E-12) / f2;

	*dp = (1.114040803E-01 - d1 * 7.739E-04 * 2.0E-08 - d12 * 2.6E-07 * 3.0E-12) / f2;

	*dnp = (5.29539222E-02 - d1 * 1.557E-04 * 2.0E-08 - d12 * 5.0E-08 * 3.0E-12) / f2;
}

/*****************************************************************************/
int basis2d (struct WebTideInfo *W, double ptx, double pty, double *basis, int *in)
/* Finds the W->closest node to a point (ptx, pty) and the element containing
   that point, if one exists.
   Also gets the basis functions for interpolations to that point.
   Returns the element number containing the point.
   Returns -1 if no element found containing the point. */
{
	int i, n11, n22, n33;
	double xlocal[3], ylocal[3];

	/* First do a crude check on the data region */
	if (ptx < W->xmin || ptx > W->xmax || pty < W->ymin || pty > W->ymax) return (-1);

	/* Now find the W->closest node */
	W->closest = closestnode (W, ptx, pty);

	/* Try to find an element that contains the point and the W->closest node. */
	in = W->in_base + W->minelem[W->closest] * 3;
	for (i = W->minelem[W->closest]; i <= W->maxelem[W->closest]; i++) {
		n11 = *in++;
		n22 = *in++;
		n33 = *in++;
		if ((W->closest == n11) || (W->closest == n22) || (W->closest == n33)) {
			xlocal[0] = W->xdata[n11];
			ylocal[0] = W->ydata[n11];
			xlocal[1] = W->xdata[n22];
			ylocal[1] = W->ydata[n22];
			xlocal[2] = W->xdata[n33];
			ylocal[2] = W->ydata[n33];
	/* See if the point is within this element */
			if (raybound (xlocal, ylocal, ptx, pty) == 1) { /* The point is within the element */
				phi2d (xlocal, ylocal, ptx, pty, basis);
				return (i);
			}
		}
	}

#ifdef SEARCH_ALL
	/* If the W->closest node's elements don't work, search through all elements */
	in = W->in_base;
	for (i = 0; i < W->numelems; i++) {
		n11 = *in++;
		n22 = *in++;
		n33 = *in++;
		xlocal[0] = W->xdata[n11];
		ylocal[0] = W->ydata[n11];
		xlocal[1] = W->xdata[n22];
		ylocal[1] = W->ydata[n22];
		xlocal[2] = W->xdata[n33];
		ylocal[2] = W->ydata[n33];
	/* See if the point is within this element */
		if (raybound (xlocal, ylocal, ptx, pty) == 1) {
			phi2d (xlocal, ylocal, ptx, pty, basis);
			return (i);
		}
	}
#endif

	return (-1);
}

/*****************************************************************************/
char *caseup (char *lower)
/* Convert a string to all upper case. */
{
	char  toup[1], unique[50] = {0};
	int leng, i;

	leng = strlen (lower);
	for (i=0; i<leng; i++) {
		sscanf (lower, "%c", toup);;
		lower++;
		unique[i] = toupper (toup[0]) ;
	}
	unique[leng] ='\0';

	return (strdup (unique));
}

/*****************************************************************************/
int closestnode (struct WebTideInfo *W, double ptx, double pty)
/* Find the node that is closest to the point (ptx, pty). */
{
	int close = 0, i;
	double currdist, closedist = 1e30, coslat;

	coslat = cos (pty * DEG2RAD);
	for (i = 0; i < W->numnodes; i++) {
/* Change to Manhattan distance - Aug 7 2003*/
		currdist = dis_man (W->xdata[i], W->ydata[i], ptx, pty, coslat);
/*  currdist = dis_sq (W->xdata[i], W->ydata[i], ptx, pty, coslat); */
		if (currdist < closedist) {
			closedist = currdist;
			close = i;
		}
	}
	return (close);
}

/*****************************************************************************/
double dis_sq (double lng1, double lat1, double lng2, double lat2, double coslat)
/* Calculate the distance between to latitude-longitude points */
{
	double x, y, r;

	x = coslat * (lng1 - lng2);
	y = lat1 - lat2;
	r = x * x + y * y;
	return (r);
}

/*****************************************************************************/
double dis_man (double lng1, double lat1, double lng2, double lat2, double coslat)
/* Calculate the "Manhattan" distance between to latitude-longitude points */
{
	double r;

	r = fabs (lng1 - lng2) + fabs (lat1 - lat2);
	return (r);
}

/*****************************************************************************/
int openvuf (struct WebTideInfo *W, char *iosfile)
/* Read in the constituent data from the IOS_tidetbl file */
{
	FILE *VUF;
	int cnt, d1, d2, d3, d4, d5, d6, i, j, nln, nsat, num;
	double fact1, fact2, fact3, fact4, phse;
	char *dummy, *inp, *name1, *name2, *name3, *name4;
	char *satin1, *satin2, *satin3;
	struct mainptr *newcon;
	struct shallptr *newshall;
	struct sconptr *newscon;

/* Allocate memory for input variables */
	inp = (char *)malloc (92 * sizeof (char));
	dummy = (char *)malloc (21 * sizeof (char));
	name1 = (char *)malloc (21 * sizeof (char));
	name2 = (char *)malloc (21 * sizeof (char));
	name3 = (char *)malloc (21 * sizeof (char));
	name4 = (char *)malloc (21 * sizeof (char));
	satin1 = (char *)malloc (25 * sizeof (char));
	satin2 = (char *)malloc (25 * sizeof (char));
	satin3 = (char *)malloc (25 * sizeof (char));

	if ( (VUF = fopen (iosfile, "r")) == NULL) {
		fprintf (stderr, "Could not open >IOS_tidetbl<!!!\n");
		fprintf (stderr, "Exiting ... \n\n");
		exit (-11);
	}

/* Counters for the number of main and shallow water constituents */
	W->nmain = 0;
	W->nshall = 0;

/* Read in the main constituents*/
	while (fgets (inp, 90, VUF) != NULL) {
		if (strlen (inp) < (unsigned)name_len) {
			break;
/* A blank line denotes the end of the main constituents */
		}

		sscanf (inp, "%s %d %d %d %d %d %d %lf %d", name1, &d1, &d2, &d3, &d4, &d5, &d6, &phse, &nsat);
		W->nmain++;

/* Re-allocate memory if necessary */
		if (W->nmain > maxcons) {
			W->cons = (struct mainptr**)realloc (W->cons, W->nmain * sizeof (struct mainptr*));
		}
		if (W->cons == NULL) {
			fprintf (stderr," Error in allocating memory for main constituents.\n");
			fprintf (stderr,  "Exiting ... \n\n");
			exit (-21);
		}

/* Create a new main constituent node */
		newcon = (struct mainptr *)malloc (sizeof (struct mainptr));
		if (newcon) {
			newcon->name = (char *)malloc (name_len * sizeof (char));
			strncpy (newcon->name, name1, name_len);
			newcon->dood[0] = d1;
			newcon->dood[1] = d2;
			newcon->dood[2] = d3;
			newcon->dood[3] = d4;
			newcon->dood[4] = d5;
			newcon->dood[5] = d6;
			newcon->phase = phse;
			newcon->nsats = nsat;

/* Read in the satellites for this constituent, if any */
			if (nsat > 0) {
				j = 0;
				newcon->sats = (struct satptr **)malloc (nsat * sizeof (struct satptr *));
				nln = ( (nsat - 1) / 3) + 1;
				for (i = 0; i < nln; i++) {
					if (fgets (inp, 90, VUF) == NULL) {
						fprintf (stderr,  "Error in reading %s satellite.\n", newcon->name);
						fprintf (stderr,  "Exiting ... \n\n");
						exit (-12);
					}
					cnt = nsat - (i * 3); /* # of satellites on this line */
					switch (cnt) {
						case 1:
							sscanf (inp, "%12c%23c", dummy, satin1);
							newcon->sats[j] = add_sat (satin1);
							j++;
							break;
						case 2:
							sscanf (inp, "%12c%23c%23c", dummy, satin1, satin2);
							newcon->sats[j] = add_sat (satin1);
							j++;
							newcon->sats[j] = add_sat (satin2);
							j++;
							break;
						default:
							sscanf (inp, "%12c%23c%23c%23c", dummy, satin1, satin2, satin3);
							newcon->sats[j] = add_sat (satin1);
							j++;
							newcon->sats[j] = add_sat (satin2);
							j++;
							newcon->sats[j] = add_sat (satin3);
							j++;
							break;
					}
				}
			}
			else {
				/*  NO satellites */
				newcon->sats = NULL;
			}
			W->cons[W->nmain-1] = newcon;
		}
	}

/* Read in the shallow water constiuents */
	while (fgets (inp, 90, VUF) != NULL) {
/* A blank line denotes the end of the shallow water constituents */
		if (strlen (inp) < (unsigned)name_len) break;

		sscanf (inp, "%s %d", name1, &num);
		(W->nshall)++;
/* Re-allocate memory if necessary */
		if (W->nshall > maxshll) W->shall = (struct shallptr**)realloc (W->shall, W->nshall * sizeof (struct shallptr*));
		if (W->shall == NULL) {
			fprintf (stderr," Error in allocating memory for shallow water const.\n");
			fprintf (stderr, "Exiting ... \n\n");
			exit (-21);
		}

/* Create a new shallow water constituent node */
		newshall = (struct shallptr *)malloc (sizeof (struct shallptr));
		if (newshall) {
			newshall->name = (char *)malloc (name_len * sizeof (char));
			strncpy (newshall->name, name1, name_len);
			strcpy (name1, ""); strcpy (name2, "");
			strcpy (name3, ""); strcpy (name4, "");
			newshall->numcon = num;
			newshall->shcon = (struct sconptr **)malloc (num * sizeof (struct sconptr *));
			switch (num) {  /* # of main constituent factors */
/* For each factor, create a node, get the info and add it to the shallow water constituent node */
				case 4:
					sscanf (inp, "%s %d %lf %s %lf %s %lf %s %lf %s", dummy, &i, &fact1, name1, &fact2, name2, &fact3, name3, &fact4, name4);
	  				newscon = (struct sconptr *)malloc (sizeof (struct sconptr));
					newscon->name = (char *)malloc (name_len * sizeof (char));
	  				strncpy (newscon->name, name4, name_len);
	  				newscon->factor = fact4;
	  				newshall->shcon[3] = newscon;
				case 3:
					sscanf (inp, "%s %d %lf %s %lf %s %lf %s", dummy, &i, &fact1, name1, &fact2, name2, &fact3, name3);
	  				newscon = (struct sconptr *)malloc (sizeof (struct sconptr));
					newscon->name = (char *)malloc (name_len * sizeof (char));
	  				strncpy (newscon->name, name3, name_len);
	  				newscon->factor = fact3;
	  				newshall->shcon[2] = newscon;
				case 2:
					if (strlen (name2) == 0) sscanf (inp, "%s %d %lf %s %lf %s", dummy, &i, &fact1, name1, &fact2, name2);
	  				newscon = (struct sconptr *)malloc (sizeof (struct sconptr));
					newscon->name = (char *)malloc (name_len * sizeof (char));
	  				strncpy (newscon->name, name2, name_len);
	  				newscon->factor = fact2;
	  				newshall->shcon[1] = newscon;
				case 1:
					if (strlen (name1) == 0) sscanf (inp, "%s %d %lf %s", dummy, &i, &fact1, name1);
	  				newscon = (struct sconptr *)malloc (sizeof (struct sconptr));
					newscon->name = (char *)malloc (name_len * sizeof (char));
	  				strncpy (newscon->name, name1, name_len);
	  				newscon->factor = fact1;
	  				newshall->shcon[0] = newscon;
					break;
			}
			W->shall[W->nshall-1] = newshall;
		}
	}

	fclose (VUF);

/* Free the memory of the input variables */
	free (inp);
	free (dummy);
	free (name1);
	free (name2);
	free (name3);
	free (name4);
	free (satin1);
	free (satin2);
	free (satin3);

	return (0);
}

/*****************************************************************************/
void phi2d (double *xloc, double *yloc, double ptx, double pty, double *basis)
/* Calculates the basis functions for interpolating to a point inside an element. */
{
	double area, a, b, c;
	int i, j, k;

	area = 0.5 * (xloc[0] * (yloc[1] - yloc[2]) + xloc[1] * (yloc[2] - yloc[0]) + xloc[2] * (yloc[0] - yloc[1]));

	/* Calculate the Basis function... */
	for (i = 0; i <= 2; i++) {
		switch (i) {
			case 0 : j = 1;  k = 2;  break;
			case 1 : j = 2;  k = 0;  break;
			case 2 : j = 0;  k = 1;  break;
		}
		a = (xloc[j] * yloc[k] - xloc[k] * yloc[j]) / (area * 2);
		b = (yloc[j] - yloc[k]) / (area * 2);
		c = - (xloc[j] - xloc[k]) / (area * 2);
		basis[i] = a + b * ptx + c * pty;
	}
}

/******************************************************************************/
int raybound (double *xd, double *yd, double ptx, double pty)
/*  Subroutine to check wether or not a point is inside a polygon.
		The process is as follows:
		- Use an arbitrary ray (here, y = constant and x >= xref), starting from
			the point and going off to infinity.
		- Count the number of polygon boundaries it crosses.
		- If an odd number, the point is inside the polygon, otherwise it is outside.
		Returns 1 if inside, 0 if outside.
*/
{
	int i, j, bcross;
	double b, m, x;

	bcross = 0; /* Number of boundary crossings. */

/* Check to see if the element side crosses the International Dateline
   (changes sign at +180/-180 degrees) and if so, change the longitudes
   so that they all have the same sign. */
  
	if( ptx > 0.0 ) { 
		if(( xd[0] < -170.0 ) && ((xd[1] > 170.0 ) || ( xd[2] > 170.0 ))) xd[0] += 360.0;
		if(( xd[1] < -170.0 ) && ((xd[0] > 170.0 ) || ( xd[2] > 170.0 ))) xd[1] += 360.0;
		if(( xd[2] < -170.0 ) && ((xd[1] > 170.0 ) || ( xd[0] > 170.0 ))) xd[2] += 360.0;
	} else {
		if(( xd[0] > 170.0 ) && ((xd[1] < -170.0 ) || ( xd[2] < -170.0 ))) xd[0] -= 360.0;
		if(( xd[1] > 170.0 ) && ((xd[0] < -170.0 ) || ( xd[2] < -170.0 ))) xd[1] -= 360.0;
		if(( xd[2] > 170.0 ) && ((xd[1] < -170.0 ) || ( xd[0] < -170.0 ))) xd[2] -= 360.0;
	}
  
/* As above, except for the Greenwich meridian, for longitude coordinates
   that range from 0 to 360 degrees. */

	if( ptx > 350.0 ) {
		if(( xd[0] < 10.0 ) && ((xd[1] > 350.0 ) || ( xd[2] > 350.0 ))) xd[0] += 360.0;
		if(( xd[1] < 10.0 ) && ((xd[0] > 350.0 ) || ( xd[2] > 350.0 ))) xd[1] += 360.0;
		if(( xd[2] < 10.0 ) && ((xd[1] > 350.0 ) || ( xd[0] > 350.0 ))) xd[2] += 360.0;
	} else {
		if(( xd[0] > 350.0 ) && ((xd[1] < 10.0 ) || ( xd[2] < 10.0 ))) xd[0] -= 360.0;
		if(( xd[1] > 350.0 ) && ((xd[0] < 10.0 ) || ( xd[2] < 10.0 ))) xd[1] -= 360.0;
		if(( xd[2] > 350.0 ) && ((xd[1] < 10.0 ) || ( xd[0] < 10.0 ))) xd[2] -= 360.0;
	}

	for (i = 0; i <= 2; i++) {

	/* for each line segment around the element */
		j = ( (i == 2) ? 0 : i + 1);

	/* If both endpoints of the line segment are on the same (vertical)
		side of the ray, do nothing.
		Otherwise, count the number of times the ray intersects the segment. */

		if (!(( (yd[i] < pty) && (yd[j] < pty)) ||
				 ( (yd[i] >= pty) && (yd[j] >= pty)))) {

			if (xd[i] != xd[j]) {
				m = (yd[j] - yd[i]) / (xd[j] - xd[i]);
				b = yd[i] - m * xd[i] ;
				x = (pty - b) / m ;
				if (x > ptx) bcross++;
			}
			else if (xd[j] > ptx)
				bcross++;
		}
	}

/*  Return the evenness/oddness of the boundary crossings
				i.e. the remainder from division by two. */
	return (bcross % 2);
}

/*****************************************************************************/
void readnodes (struct WebTideInfo *W, char *nfile, char *efile)
/* Read in the mesh from a .nod/.ele set of files. */
{
	FILE *fdn, *fde;
	double xdum, ydum;
	int i, *in, in1, in2, in3, dum, zero = 0;

	/* First get the number of nodes in the mesh. */

	fdn = fopen (nfile, "r");
	W->numnodes = -1;
	while (!feof (fdn)) {
		W->numnodes++;
		fscanf (fdn, "%d %lf %lf", &dum, &xdum, &ydum);
	}

	/* Get the number of elements in the mesh. */

	fde = fopen (efile, "r");
	W->numelems = -1;
	while (!feof (fde)) {
		W->numelems++;
		fscanf (fde, "%d %d %d %d", &dum, &in1, &in2, &in3);
	}

	/* Read in the nodes. */

	rewind (fdn);
	W->xdata = (double *)malloc ((W->numnodes) * sizeof (double));
	W->ydata = (double *)malloc ((W->numnodes) * sizeof (double));
	W->minelem = (int *)malloc ((W->numnodes) * sizeof (int));
	W->maxelem = (int *)malloc ((W->numnodes) * sizeof (int));
	W->xmin = 1e30, W->xmax = -1e30, W->ymin = 1e30, W->ymax = -1e30;
	for (i = 0; i < W->numnodes; i++) {
	/* The program starts counting at zero. Does the data file? */
		fscanf (fdn, "%d %lf %lf", &dum, &W->xdata[i], &W->ydata[i]);
		if (i == 0) zero = dum;
		if (W->xdata[i] < W->xmin) W->xmin = W->xdata[i];
		if (W->xdata[i] > W->xmax) W->xmax = W->xdata[i];
		if (W->ydata[i] < W->ymin) W->ymin = W->ydata[i];
		if (W->ydata[i] > W->ymax) W->ymax = W->ydata[i];
		W->minelem[i] = W->numelems;
		W->maxelem[i] = 0;
	}
	fclose (fdn);

	/* - Read in the elements.
		 - Correct for the starting number in the nodal data set (zero)
		 - Tally the first and last element containing each node
	*/

	rewind (fde);
	W->in_base = (int *)malloc (3 * W->numelems * sizeof (int));
	in = W->in_base;
	for (i = 0; i < W->numelems; i++) {
		fscanf (fde, "%d %d %d %d", &dum, &in1, &in2, &in3);
		in1 -= zero; in2 -= zero; in3 -= zero;
		if (i < W->minelem[in1]) W->minelem[in1] = i;
		if (i > W->maxelem[in1]) W->maxelem[in1] = i;
		if (i < W->minelem[in2]) W->minelem[in2] = i;
		if (i > W->maxelem[in2]) W->maxelem[in2] = i;
		if (i < W->minelem[in3]) W->minelem[in3] = i;
		if (i > W->maxelem[in3]) W->maxelem[in3] = i;
		*in++ = in1; *in++ = in2; *in++ = in3;
	}
	fclose (fde);
}

/*****************************************************************************/
void setvuf (struct WebTideInfo *W, long int kd, double xlat)
/* Calculate the amplitudes, phases, etc. for each of the constituents */
{
	double d1, tau, dtau, slat, sumc, sums, v, vdbl;
	double adj = 0, dum, dd[6], uu, uudbl;
	double h, pp, s, p, np, dh, dpp, ds, dp, dnp;
	int i, j, k;
	struct satptr *sat;

	d1 = (double)(kd) - 0.5;

	astro_angles (d1, &h, &pp, &s, &p, &np, &dh, &dpp, &ds, &dp, &dnp);

	tau = h - s;
	dtau = 365.00 + dh - ds;
	slat = sin (xlat * DEG2RAD);

/* The main constituents */
	for (k = 0; k < W->nmain; k++) {
		for (i = 0; i < 6; i++) dd[i] = W->cons[k]->dood[i];
		W->cons[k]->freq = (dd[0] * dtau + dd[1] * ds + dd[2] * dh + dd[3] *dp + dd[4] * dnp + dd[5] * dpp) / (24.0 * 365.0);
		vdbl = dd[0] * tau + dd[1] * s + dd[2] * h + dd[3] * p + dd[4] * np + dd[5] * pp + W->cons[k]->phase;
/*    v = vdbl - (double) ((int)(vdbl) / 2 * 2.0);*/
		v = vdbl - (floor (floor (vdbl) / 2.0) * 2.0);

		sumc = 1.0;
		sums = 0.0;

		for (i = 0; i < W->cons[k]->nsats; i++) {
			sat = W->cons[k]->sats[i];
			switch (sat->corr) {
				case 0: adj = sat->ratio; break;
				case 1: adj = sat->ratio * 0.36309 * (1.0 - 5.0 * slat * slat) / slat; break;
				case 2: adj = sat->ratio * 2.59808 * slat; break;
			}
			uudbl = (double)sat->deld[0] * p + (double)sat->deld[1] * np + (double)sat->deld[2] * pp + (double)sat->phase;
			uu = modf (uudbl, &dum);
			sumc += (adj * cos (uu * TWO_PI));
			sums += (adj * sin (uu * TWO_PI));
		}
		W->cons[k]->f = sqrt ( (sumc * sumc) + (sums * sums));
		W->cons[k]->vu = v + atan2 (sums, sumc) / TWO_PI;
	}

/* The shallow water constituents */
	for (k = 0; k < W->nshall; k++) {
		W->shall[k]->f = 1.0;
		W->shall[k]->vu = 0.0;
		W->shall[k]->freq = 0.0;
		for (i = 0; i < W->shall[k]->numcon; i++) {
			for (j = 0; j < W->nmain; j++) {
				if (strcmp (W->cons[j]->name, W->shall[k]->shcon[i]->name) == 0) {
					W->shall[k]->f *= pow (W->cons[j]->f, fabs (W->shall[k]->shcon[i]->factor));
	  W->shall[k]->vu += (W->shall[k]->shcon[i]->factor * W->cons[j]->vu);
	  W->shall[k]->freq += (W->shall[k]->shcon[i]->factor * W->cons[j]->freq);
	  break;
	}
			}
		}
	}
}

/*****************************************************************************/
double TideP (struct WebTideInfo *W, long int kd, double dthr, double latitude, double *ampl, double *phase)
/* Calculates and returns the tidal correction */
{
	int i, indx;
	double dum, radgmt, revgmt, res;

/*  kd--;*/
	setvuf (W, kd, latitude);
	res = 0.0;

/* For each of the desired constituents ... (See top of program) */
	for (i = 0; i < W->nconst; i ++) {
/* Find the constituent from those loaded from IOS_tidetbl */
		indx = vuf (W, W->constnames[i]);
		if (indx < 0) {
			fprintf (stderr, "Bad Input Constituent: %s\n", W->constnames[i]);
			fprintf (stderr, "Exiting ... \n\n");
			exit (-1);
		}

		if (indx < W->nmain) {                        /* Main constituent */
			revgmt = W->cons[indx]->freq * dthr + W->cons[indx]->vu - phase[i] / 360.0;
			radgmt = TWO_PI * modf (revgmt, &dum);
			res += W->cons[indx]->f * ampl[i] * cos (radgmt);
		}
		else if ( (indx - W->nmain) < W->nshall) {     /* Shallow water constituent */
			indx -= W->nmain;
			revgmt = W->shall[indx]->freq * dthr + W->shall[indx]->vu - phase[i] / 360.0;
			radgmt = TWO_PI * modf (revgmt, &dum);
			res += W->shall[indx]->f * ampl[i] * cos (radgmt);

		}
		else {
			fprintf (stderr, "Error in index: %d %d %d\n", indx, W->nmain, W->nshall);
			fprintf (stderr, "Exiting ... \n\n");
			exit (-2);
		}
	}
	return (res);
}

/*****************************************************************************/
int vuf (struct WebTideInfo *W, char *inname)
/* Finds constituent info corresponding to inname and returns the index to
   the node containing the info. Shallow water constituent indices are
   returned as their number greater than the max # of main constituents.
   e.g. if we want the 2nd shallow water constituent and there are 45
   main constituents, then the index returned is 46, since constituents
   are counted from zero. (45 - 1 + 2 = 46) */
{
	int i, j = 0;

	i = 0;
	while (strcmp (W->cons[i]->name, inname) != 0) {
		i++;
		if (i == W->nmain) break;
	}
	if (i == W->nmain) {
		j = 0;
		while (strcmp (W->shall[j]->name, inname) != 0) {
			j++;
			if (j == W->nshall) break;
		}
	}

	if (i < W->nmain)
		return (i);
	else if (j < W->nshall)
		return (i + j);
	else {
		fprintf (stderr, "Constituent %s not found!\n", inname);
		fprintf (stderr, "Exiting ... \n\n");
		exit (-9);
	}
}
