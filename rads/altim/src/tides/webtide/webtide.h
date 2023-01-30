/* Header file for webtide subroutines */

/* The following are structures used to hold the info of the constituents
   First is the structure to hold the satellites for each main constituent */
struct satptr {
	int deld[3];
	double phase;
	double ratio;
	int corr;
};

/* This is for the main constituents */
struct mainptr {
	char *name;
	int dood[7];
	double phase;
	int nsats;
	struct satptr **sats;
	double f;
	double freq;
	double vu;
};

/* This holds the main constituent factors for the shallow water constituents */
struct sconptr {
	char *name;
	double factor;
};

/* This structure holds the info for the shallow water constituents */
struct shallptr {
	char *name;
	int numcon;
	struct sconptr **shcon;
	double f;
	double freq;
	double vu;
};

/* Structure to hold all the global variables */
struct WebTideInfo
{
	long int nconst;	/* # of constituents */
	double *xdata, *ydata;	/* Pointers to longitude and latitude */
	int *minelem, *maxelem;	/* Lowest and highest element containing a certain node */
	double xmin, xmax, ymin, ymax;	/* Model boundaries */
	int closest;		/* Closest node number to an arbitrary point */
	int *in_base;	/* Pointer to the input data */
	int numnodes, numelems;	/* Total # of nodes and elements */
	char **constnames;
	int cplx;		/* Is input data complex (velocities)? */
	struct mainptr **cons;
	struct shallptr **shall;
	double *amp_base, *phase_base;	/* Pointers to amplitudes and phases */
	double *amp2_base, *phase2_base;	/* Idem, 2nd component */
#ifdef MIN_TIDE
	double *min_tide_base, *min_tide2_base; /* Minimum tide values */
#endif
	int nmain, nshall;
};

long int WebTideInit (char *dirname, char *ext, struct WebTideInfo *W);
long int WebTide (struct WebTideInfo *W, double *time, double *latitude, double *longitude, double *tide, double *tide2);
void WebTideFree (struct WebTideInfo *W);
