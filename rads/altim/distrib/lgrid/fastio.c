/*FASTIO -- Fast low-level UNIX I/O routines
*+
*
* FASTIO is a set of functions that makes fast, low-level Unix I/O routines
* available to a Fortran program, accounting for some of the differences in the
* C/Fortran interfaces of different machines.
*
* Specifically, some linkers expect a C-routine which is called from a 
* Fortran program to have a name ending in an underscore.  So, for example,
* if the Fortran program calls "openf()", the C-function's name must
* be "openf_()".  Other vendor's (like HP and RS6000) look for a C-routine with
* the actual name as called from Fortran (without the underscore).
* 
* Therefore, one can define UNDERSCOREAFTER during the compilation to get
* the underscores after the function names. Obviously, one can compile
* with and without the -DUNDERSCOREAFTER option and link both objects
* if you are not sure. Other options are: -DUNDERSCOREBEFORE and -DCAPITALS.
*
* Secondly, when character strings are passed from Fortran to C, the
* string's length is implicitly passed, unbeknownst to the Fortran caller.
* Some C/Fortran interfaces put the string length as the next argument after
* the string itself; others place the string length at the end of the
* argument list.  Again, to support both implementations, when a string is
* passed to the following routines it is the last argument in the Fortran
* call, so that the string and its length are the last two arguments in the
* corresponding C-function.
*
* For the definition of the open flag in the call "openf()", it is required
* to include the file "fcntl.h". On most operating systems, this file
* is in "/usr/include", but we have encountered systems where it is in
* "/usr/include/sys". To support both, you should compile this source-code
* with the additional search path specified: -I/usr/include/sys.
*
* The syntax of the Fortran calls are provided below.
*
* -------
* $Log: fastio.c,v $
* Revision 1.12  2006/10/10 17:36:54  rads
* - Introduced new OPENF function
*
* Revision 2.3  2005/05/16 00:21:27  rads
* - Added include string.h
*
* Revision 2.2  2004/08/28 15:29:01  remko
* Reintroduced unistd.h (required by seekf on Mac OS X); spiffed up manual
*
* 24-Jul-2003 - Added IOFLAG function and removed O_TRUNC from read/write flag.
* 19-Jun-2003 - Included understanding of '-' to indicate stdin/stdout.
* 20-Dec-2001 - Improved manual, added O_TRUNC to the read/write flag.
* 12-Oct-2001 - Added line reading and writing. Add read/write flag to IOCONST.
*  7-Jan-1996 - Underscores through UNDERSCOREAFTER. Warning also without
*               carriage return.
* 17-May-1994 - Nice manual
* 11-Nov-1993 - Addition of seekf and warning by Remko Scharroo, DUT/SSR&T
*  2-Dec-1992 - fastio.c: John L. Lillibridge, NOAA/NOS/OES Geosciences Lab
* -------

*&OPENF -- Open file for FASTIO
*+
*     FUNCTION OPENF (FNAME, MODE)
*     CHARACTER*(*) FNAME, MODE
*
* OPENF opens file FNAME for reading and or writing, depending on the
* string MODE, which is the same as the "mode" argument of the C function
* fopen().
*
* FNAME can be up to 1024 characters long. Specify '-' for reading from
* standard input or writing to standard output.
*
* MODE is a 1- or 2-byte string indicating whether the file needs to
* be opened for reading, writing, or both, and whether trunctation needs
* to happen. The different options are:
* 'r'  : Open file for reading; start at the beginning of the file.
* 'r+' : Open file for reading and writing; start at the beginning
*        of the file.
* 'w'  : Truncate file to zero length or create it for writing only.
* 'w+' : Open file for reading and writing. Truncate it to zero length
*        if it exists, create it otherwise.
* 'a'  : Open for writing. Create the file if it does not exist. Start
*        at the end of the file. Subsequent writes are appended at the end.
* 'a+' : Open for reading and writing. Create the file if it does not
*        exist. Start at the end of the file. Subsequent write are
*        appended at the end of the file.
*
* OPENF returns the file descriptor for use in subsequent calls to
* readf, writef, or closef. If OPENF is negative, an error occurred while
* opening the file. Use PERRORF to report the error.
*
* Usage for reading an EXISTING file.
*
*     INTEGER FD, OPENF, STRING
*
*     FD = OPENF ('input_file', 'r')
*     CALL READF (FD, 4, STRING)
*     CALL CLOSEF (FD)
*
* Usage for writing a file which may or may not exist.
*
*     FD = OPENF ('output_file', 'w')
*     CALL WRITEF (FD, 4, STRING)
*     CALL CLOSEF (FD)
*
* Finally, to read and write to an existing file, use:
*
*     FD = OPENF ('io_file', 'r+')
*
* Note: Standard input and standard output can be opened through OPENF
* by specifying '-' as the file name FNAME and using MODE 'r' and 'w'
* respectively
*
* Input argument:
*  FNAME : File name of the input or output file. Maximum 1024 chars.
*          Use '-' for standard input or standard output.
*  MODE  : Open mode for reading or writing, see above.
*
* Returned value:
*  OPENF : Contains the file descriptor on return.
*          If the returned value is negative an error occurred while
*          opening the file.
*-

*&CLOSEF -- Close file from FASTIO access
*+
*     FUNCTION CLOSEF (FD)
*     INTEGER CLOSEF (FD)
*
* Closes the file with descriptor FD from FASTIO access. CLOSEF returns
* 0 when properly closed. Otherwise, use PERRORF to report the error.
* 
* Usage:
*      IOS = CLOSEF (FD)
* or:
*      CALL CLOSEF (FD)
*
* In the last case the return code is ignored.
*
* Input argument:
*  FD     : File descriptor returned by OPENF.
* Returned value:
*  CLOSEF : Error code or 0 on proper closing.
*-

*&READF -- FASTIO read routine
*+
*      FUNCTION READF (FD, NBYTE, BUFFER)
*      INTEGER FD, NBYTE, READF
*      BYTE    BUFFER (NBYTE)
*
* Reads a string of bytes or a line (until carriage return) from the file
* associated with descriptor FD (returned by the OPENF call).
* When reading from standard input, use FD=0.
*
* When reading a line up to AND INCLUDING the carriage return, use NBTYE=0.
* In this case the string that is read may contain no more than 1024 bytes
* and will end with a line-feed (CHAR(10)). Upon return, the value of
* READF will be the number of bytes of the string, including the line-feed.
*
* When reading a string of a specified length (binary or otherwise), NBYTE
* specifies the length of the string to be read. READF will reflect the
* number of bytes actually read.
*
* The array BUFFER will hold the data, but can (of course) also be
* associated with any other string, scalar, or n-dimensional array.
*
* The function returns the number of bytes actually read in READF. If
* READF < 0, a read error occurred which can be reported by calling
* PERRORF. READF = 0 indicates an end-of-file mark.
*
* Example of binary reading:
*
*      INTEGER FD, RBYTE, READF, BUFFER(20)
*      RBYTE = READF (FD, 80, BUFFER)
*      IF (RBYTE .EQ. 80) THEN
*         WRITE (6,*) 'All bytes returned'
*         WRITE (6,*) BUFFER
*      ELSE IF (RBYTE .EQ. 0) THEN
*         WRITE (6,*) 'End-of-file'
*      ELSE IF (RBYTE .LT. 0) THEN
*         CALL PERRORF ('readf error')
*      ELSE
*         WRITE (6,*) 'Premature EOF after ',RBYTE,' bytes'
*      ENDIF
*
* Example of line reading:
*
*      INTEGER FD, RBYTE, READF
*      CHARACTER*80 LINE
*      RBYTE = READF (FD, 0, LINE)
*      IF (RBYTE .EQ. 0) THEN
*         WRITE (6,*) 'End-of-file'
*      ELSE IF (RBYTE .LT. 0) THEN
*         CALL PERRORF ('readf error')
*      ELSE
*         WRITE (6,*) 'Read ',RBYTE,' bytes'
*         WRITE (6,*) 'String = ',LINE(:RBYTE-1)   ! Strip off trailing CR
*      ENDIF
*
* Input arguments:
*  FD     : File descriptor returned by OPENF.
*  NBYTE  : Number of bytes to be read. Use NBYTE=0 to read to next carriage return.
*
* Output argument:
*  BUFFER : Buffer containing the bytes that have been read.
*
* Returned value:
*  READF  : Number of bytes read, or (if negative) error code.
*-

*&WRITEF -- FASTIO write routine
*+
*     FUNCTION WRITEF (FD, NBYTE, BUFFER)
*     INTEGER FD, NBYTE, WRITEF
*     BYTE    BUFFER(NBYTE)
*
* Writes a string of bytes or a line (followed by carriage return) into
* the file associated by descriptor FD (which is returned by the OPENF
* call). When writing to standard output use FD=1.
*
* When a line is written to the file, simply specify NBYTE=0. The string
* should then include a trailing line-feed (CHAR(10)) and may not be
* longer that 1024 bytes.
* The function value of WRITEF will be the number of bytes of the string
* including the line-feed.
*
* When writing a binary string, NBYTE specifies the number of bytes
* contained of the binary string. WRITEF will reflect the actual
* numbers of bytes written.
*
* The array BUFFER contains the data that has to be written, but can
* (of course) also be associated with any other string, scalar, or
* n-dimensional array.
*
* The function returns the number of bytes actually written in WRITEF. If
* WRITEF < 0, a write error occurred which can be reported by calling
* PERRORF.
*
* Input arguments:
*  FD     : File descriptor returned by OPENF
*  NBYTE  : Number of bytes to be written. Use 0 for a string ending in
*           carriage return.
*  BUFFER : Buffer containing the bytes that have to be written
*
* Returned value:
*  WRITEF : Number of bytes written, or (if negative) error code.
*-

*&SEEKF -- FASTIO file positioning routine
*+
*     FUNCTION SEEKF (FD, OFFSET, WHENCE)
*     INTEGER FD, OFFSET, WHENCE, SEEKF
*
* Spools file to the position indicated by the parameter OFFSET (in
* bytes). The parameter WHENCE indicates whether OFFSET is counted
* - from the beginning of the file (WHENCE = 0),
* - from the current position      (WHENCE = 1),
* - from the end of the file       (WHENCE = 2)
* Upon return SEEKF reports the number of bytes skipped, or contains
* an error code (to be reported by PERRORF).
*
* Input arguments:
*  FD     : File descriptor returned by OPENF
*  OFFSET : Offset in bytes from the point indicated by WHENCE
*  WHENCE : Indicates whether the OFFSET is counted from the
*           start of the file (0), current position (1), or the
*           end of the file (2).
*
* Returned value:
*  SEEKF  : On proper return: new position in file, otherwise error code
*-

*&PERRORF -- Report last FASTIO error
*+
*     SUBROUTINE PERRORF (COMMENT)
*     CHARACTER*(*) COMMENT
*
* Reports the last low-level I/O error returned by the system. The form of
* the reported string which is sent to standard-error is:
* "COMMENT: error message", where COMMENT is held in the parameter in the
* call to PERRORF.
*
* Input argument:
*  COMMENT : Additional comment to the error message
*-
*=
*/

#define stringlength 1024

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>

#include <sysdep.h>
#ifdef CAPITALS
#define openf OPENF
#define closef CLOSEF
#define readf READF
#define writef WRITEF
#define perrorf PERRORF
#define seekf SEEKF
#endif

#ifdef UNDERSCOREBEFORE
#define openf _openf
#define closef _closef
#define readf _readf
#define writef _writef
#define perrorf _perrorf
#define seekf _seekf
#endif

#ifdef UNDERSCOREAFTER
#define openf openf_
#define closef closef_
#define readf readf_
#define writef writef_
#define perrorf perrorf_
#define seekf seekf_
#endif
/*****************************************************************************/

long openf(char *fname, char *mode, int slen, int mlen)
{
	char temp[stringlength];	/* Scratch string...could malloc() */
	int len = slen - 1;
	strncpy(temp,fname,slen);

/* Strip off trailing blanks and null terminate string */

	while (temp[len] == ' ' && len >= 0) len--;
	temp[len+1] = '\0';
 
	if (len == 0 && temp[0] == '-'){
		if (mode[0] == 'r') return (0);	/* stdin */
		else if (mode[0] == 'w' || mode[0] == 'a') return (1);	/* stdout */
	}
	else if (mlen > 1 && mode[1] == '+') {
		if (mode[0] == 'r') return(open (temp,O_RDWR,0666));
		else if (mode[0] == 'w') return(open (temp,O_RDWR|O_TRUNC|O_CREAT,0666));
		else if (mode[0] == 'a') return(open (temp,O_RDWR|O_APPEND,0666));
	}
	else if (mode[0] == 'r') return(open (temp,O_RDONLY,0444));
	else if (mode[0] == 'w') return(open (temp,O_WRONLY|O_TRUNC|O_CREAT,0666));
	else if (mode[0] == 'a') return(open (temp,O_WRONLY|O_APPEND,0666));

	fprintf(stderr, "Unknown mode in openf\n");
	return (-1);
}

/*****************************************************************************/

long closef(fd)
int *fd;
{
	if(*fd == 0) return(fclose(stdin));
	else if(*fd == 1) return(fclose(stdout));
	else return(close(*fd));
}

/*****************************************************************************/

long readf(fd,nbytes,buf,slen)
int *fd, *nbytes, slen;
void *buf;
{
	int i;
	char temp[stringlength];

	if(*fd == 0){			/* fd=0 means read from stdin */
		if(*nbytes == 0){		/* nbytes=0 signals read ascii "lines" */
			if(fgets(buf,slen,stdin)==NULL)	/* null-terminated, WITH <LF> char returned */
				return(0);
			else
				return(strlen(buf));
		}else			/* read nbytes of binary data */
			return(fread(buf,1,*nbytes,stdin));	/* Number of 1-byte items=*nbytes */
	}else{			/* read from file descriptor *fd */
		if(*nbytes == 0){		/* read ascii "lines" */
			for (i=0; i<slen; i++) {
				if (read(*fd,&temp[i],1)==0) return(0);
				if (temp[i]==10) break;
			}
			i++; temp[i]='\0';		/* null-terminated, with <LF> char returned */
			strncpy(buf,temp,i);
			return(i);
		}else			/* read nbytes of binary data */
			return(read(*fd,buf,*nbytes));
	}
}

/*****************************************************************************/

long writef(fd,nbytes,buf,slen)
int *fd, *nbytes, slen;
void *buf;
{
	int i;
	char temp[stringlength];
	if(*fd == 1){			/* fd=1 means write to stdout */
		if(*nbytes == 0){		/* nbytes=0 signals write ascii "lines" */
			strncpy(temp,buf,slen);
			for (i=0; i<slen; i++) if (temp[i]==10) break;
			i++; temp[i]='\0';
			return(fputs(temp,stdout));
		}else			/* write nbytes of binary data */
			return(fwrite(buf,1,*nbytes,stdout)); /* Number of items = *nbytes */
	}else{			/* write to file descriptor *fd */
		if(*nbytes == 0){		/* nbytes=0 signals write ascii "lines" */
			strncpy(temp,buf,slen);
			for (i=0; i<slen; i++) if (temp[i]==10) break;
			i++; temp[i]='\0';
			return(write(*fd,temp,i));
		}else			/* write nbytes of binary data */
			return(write(*fd,buf,*nbytes));
	}
}

/*****************************************************************************/

void perrorf(string,slen)
char *string;
int slen;
{
	char temp[stringlength];	/* Scratch string...could malloc() */
	int len = slen - 1;
	strncpy(temp,string,slen);

/* Strip off trailing blanks and null terminate string */

	while (temp[len] == ' ' && len >= 0) len--;
	temp[len+1] = '\0';

	perror(temp);
	return;
}

/*****************************************************************************/

long seekf(fd,nbytes,whence)
int *fd, *nbytes, *whence;
{
	if(*fd == 0) return(fseek(stdin,*nbytes,*whence));
	else return(lseek(*fd,*nbytes,*whence));
}
