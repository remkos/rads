/*
 *  File to convert ascii FES files to binary FES files
 * 
 *  File      : ascii2bin.c
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.1
 *  Date      : 6 October 2004
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#include <sys/types.h>
#endif
#include <sys/stat.h>
#include <math.h>
#include <errno.h>
#include <stdarg.h>

#define MAX_PATH	1024

#ifndef M_PI
    #define M_PI        3.14159265358979323846
#endif

typedef struct _fComplex
{
    float	re;
    float	im;
} fComplex;

static double	rad	= M_PI / 180.0;








/*
// ///////////////////////////////////////////////////////////////////////////
// This function swaps a variable
//
// Parameters:
//   a:		pointer to the variable to swap
//   b:		pointer to the variable to swap
*/
static void swap(void** a, void** b)
{
    register void *tmp = *a;

    *a = *b;
    *b = tmp;
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Verifie the system integer variables notation.
//
// Return value:
//   Returns 1 for little endian else 0 for big endian notation.
*/
static int isLittleEndian(void)
{
    static short a = 1;
    static char *p = (char *)&a;

    return(*p ? 1 : 0);
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Read a string from a stream.
//
// Parameters:
//   s:		string
//   n:		max characters to read
//   stream:	stream to read
*/
static void readStr(char *s, int n, FILE *stream)
{
    char* rc	= fgets(s, n, stream);

    if( rc == NULL )
    {
	fprintf(stderr, "Unable to read input file: %s.\n", strerror(errno));
	exit(1);
    }
}




/*
// ///////////////////////////////////////////////////////////////////////////
// Scans and formats input from a string.
//
// Parameters:
//   s:		buffer to decode
//   format:	format specifier
//   n:		number of input fields
*/
static void check(const int rc, const int n)
{
    if( rc != n )
    {
	fprintf(stderr, "Unable to decode input file.\n");
	exit(1);
    }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Read formatted data from a line.
//
// Parameters:
//   line:	line to read
//   undef:	Value of a cell nondefined in the grid.
//   y:		latitude to read
//   nX		latitude samples
//   x:		longitude to read
//   values:	The data read.
*/
static void readValues(char* line,
		       const int y,
		       const int nX,
		       int* x,
		       float* values)
{
  char	*token = strtok(line, " \t");

  while(token)
  {
    check(sscanf(token, "%f", &values[y * nX + (*x)]), 1);

    (*x)++;

    token  = strtok(NULL, " \t");
  }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Write to a file.
//
// Parameters:
//   handle:	a file handle
//   buffer
//   len:	write len bytes from buffer
*/
static void writeBuffer(int handle, void *buffer, unsigned len)
{
    if( write(handle, buffer, len) == -1 )
    {
	fprintf(stderr, "Error writting binary file: %s\n", strerror(errno));
	exit(1);
    }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// function swaps a "size" byte variable, i.e. it converts their
// representation between little and big endian (in either direction).
//
// Parameters:
//   p:		pointer to the variable to swap
//   n:		size of the variable
*/
static void swapByte(void *p, int n)
{
    register char *begin;
    register char *end;
    register char  temp;

    begin =  (char*)p;
    end   = ((char*)p) +(n - 1);
    while ( begin < end )
    {
        temp = *begin;
        *(begin++) = *end;
        *(end--)   = temp;
    }
}





/*
// ///////////////////////////////////////////////////////////////////////////
// Allocates maibn memory
//
// Parameters:
//   nitems
//   size
//
// returns a pointer to the newly allocated block. If not enough space exists
// for the new block new terminates program.
*/
static void *new(size_t nitems, size_t size)
{
  void* ptr = calloc(nitems, size);

  if( ptr == NULL )
  {
    fprintf(stderr, "Out of memory.\n");
    exit(1);
  }
  return ptr;
}







/*
// ///////////////////////////////////////////////////////////////////////////
// Main program
*/
int main(int argc, char** argv)
{
    char        buffer[MAX_PATH];
    FILE*	in;
    int		out;
    int		nX;
    int		nY;
    int		x;
    int		y		= 0;
    int		pos		= 0;
    int		tide;
    int		littleEndian	= isLittleEndian();
    int         size;
    float	latMin;
    float	lonMin;
    float	latStep;
    float	lonStep;
    float	undef;
    float	unused;
    float	defValue	= 18446744073709551616.0F;
    float*	amp;
    float*	pha;
    double	sizeOfFile;
    double	sizeOfByte;
    fComplex*	grid;

    if(argc != 4)
    {
        fprintf(stderr, "Convert ascii file to binary file.\n\n");
        fprintf(stderr, "Usage: ascii2bin ascii bin type\n\n");
        fprintf(stderr, "  ascii\tascii files to convert\n");
        fprintf(stderr, "  bin\tbinary files\n");
        fprintf(stderr, "  type\t1 if ascci file contains pure harmonic tide,\n");
	fprintf(stderr, "\telse 0 if files contains radial tide.\n");
	return 1;
    }
    
    printf("Converting %s to %s.\n", argv[1], argv[2]);
    
    if((in = fopen(argv[1], "r")) == NULL)
    {
	fprintf(stderr, "Unable to open %s: %s", argv[1], strerror(errno));
	exit(1);
    };
    
#ifdef WIN32
    if( (out = open(argv[2], O_WRONLY|O_CREAT|O_BINARY, S_IREAD|S_IWRITE)) == -1 )
#else
    if( (out = open(argv[2], O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR)) == -1 )
#endif
    {
	fprintf(stderr, "Unable to create %s: %s", argv[2], strerror(errno));
	exit(1);
    };

    check(sscanf(argv[3], "%d", &tide), 1);

    /* Reading XMIN & XMAX (X <=> Longitude) */
    readStr(buffer, MAX_PATH, in);
    check(sscanf(buffer, "%f %f", &lonMin, &unused), 2);

    /* Reading YMIN & YMAY (Y <=> Latitude) */
    readStr(buffer, MAX_PATH, in);
    check(sscanf(buffer, "%f %f", &latMin, &unused), 2);

    /* Reading DX & DY */
    readStr(buffer, MAX_PATH, in);
    check(sscanf(buffer, "%f %f", &lonStep, &latStep), 2);

    /* Reading NX & DY */
    readStr(buffer, MAX_PATH, in);
    check(sscanf(buffer, "%d %d", &nX, &nY), 2);

    /* The grids should not overlap. */
    nX--;

    /* Reading MASKX, MASKY */
    readStr(buffer, MAX_PATH, in);
    check(sscanf(buffer, "%f %f", &unused, &undef), 2);

    size = nX * nY;

    grid = new(size, sizeof(fComplex));
    amp  = new(size, sizeof(float));
    pha  = new(size, sizeof(float));

    printf("Grid properties:\n");
    printf("  xMin:\t%f\n", lonMin);
    printf("  yMin:\t%f\n", latMin);
    printf("  dX:\t%f\n", lonStep);
    printf("  dY:\t%f\n", latStep);
    printf("  nX:\t%d\n", nX);
    printf("  nY:\t%d\n", nY);

    printf("Reading ascii file...\n");

    /* For all the latitudes. */
    while(y < nY)
    {
	x = 0;

	/* For all the longitudes. */
	while(x < nX)
	{
	    /* Reading the 30 values of the amplitude. */
	    readStr(buffer, MAX_PATH, in);
	    readValues(buffer, y, nX, &x, amp);

	    x -= 30;

	    /* Reading the 30 values of the phase. */
	    readStr(buffer, MAX_PATH, in);
	    readValues(buffer, y, nX, &x, pha);
	}
	/* Reading of the redundant value for the amplitude. */
	readStr(buffer, MAX_PATH, in);

	/* Reading of the redundant value for the phase. */
	readStr(buffer, MAX_PATH, in);

	y++;
    }

    for(y = 0; y < nY; y++)
    {
	for(x = 0; x < nX; x++)
	{
	    register int ix	= y * nX + x;

	    if( amp[ix] != undef && pha[ix] != undef )
	    {
		double re = tide?
    		    amp[ix] * cos(-pha[ix] * rad):
		    amp[ix] * cos( pha[ix] * rad);
		double im = tide?
		    -amp[ix] * sin(-pha[ix] * rad):
		     amp[ix] * sin( pha[ix] * rad);
		grid[ix].re = (float) re;
		grid[ix].im = (float) im;
	    }
	    else
	    {
		grid[ix].re = defValue;
		grid[ix].im = defValue;
	    }

	    if( !littleEndian )
	    {
		swapByte(&grid[ix], sizeof(fComplex));
		swap((void*)&grid[ix].re, (void*)&grid[ix].im);
	    }
	}
    }


    printf("Writting binary file...\n");

    if( !littleEndian )
    {
	swapByte(&nX,	    sizeof(int));
	swapByte(&nY,	    sizeof(int));
	swapByte(&lonMin,   sizeof(float));
	swapByte(&latMin,   sizeof(float));
	swapByte(&lonStep,  sizeof(float));
	swapByte(&latStep,  sizeof(float));
	swapByte(&defValue, sizeof(float));
    }

    writeBuffer(out, &nX, sizeof(int));
    pos += sizeof(int);

    writeBuffer(out, &nY, sizeof(int));
    pos += sizeof(int);

    writeBuffer(out, &lonMin, sizeof(float));
    pos += sizeof(float);

    writeBuffer(out, &latMin, sizeof(float));
    pos += sizeof(float);

    /* Reading lat, lon step */
    writeBuffer(out, &lonStep, sizeof(float));
    pos += sizeof(float);

    writeBuffer(out, &latStep, sizeof(float));
    pos += sizeof(float);

    writeBuffer(out, &defValue, sizeof(float));
    pos += sizeof(float);

    printf("Size of header: %d\n", pos);

    writeBuffer(out, grid, sizeof(fComplex) * nX * nY);

    sizeOfFile = lseek(out, 0, SEEK_CUR);;
    sizeOfByte = pos + ((nY - 1) * nX + (nX - 1)) * sizeof(fComplex) + sizeof(fComplex);

    printf("Size of file: %.0f\n", sizeOfFile);

    if( sizeOfFile != sizeOfByte )
    {
	fprintf(stderr, "Incorrect file size.");
	exit(1);
    }

    fclose(in);
    close(out);

    free(grid);
    free(amp);
    free(pha);

    printf("Done.\n");

    return 0;
}


