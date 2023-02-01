/*FASTIO -- Fast low-level UNIX I/O routines
 +

  FASTIO is a set of functions that makes fast, low-level Unix I/O routines
  available to a Fortran program, accounting for some of the differences in the
  C/Fortran interfaces of different machines.

  Specifically, some linkers expect a C-routine which is called from a 
  Fortran program to have a name ending in an underscore.  So, for example,
  if the Fortran program calls "openf()", the C-function's name must
  be "openf_()".  Other vendor's (like HP and RS6000) look for a C-routine with
  the actual name as called from Fortran (without the underscore).
  
  Therefore, one can define UNDERSCOREAFTER during the compilation to get
  the underscores after the function names. Obviously, one can compile
  with and without the -DUNDERSCOREAFTER option and link both objects
  if you are not sure. Other options are: -DUNDERSCOREBEFORE and -DCAPITALS.

  Secondly, when character strings are passed from Fortran to C, the
  string's length is implicitly passed, unbeknownst to the Fortran caller.
  Some C/Fortran interfaces put the string length as the next argument after
  the string itself; others place the string length at the end of the
  argument list.  Again, to support both implementations, when a string is
  passed to the following routines it is the last argument in the Fortran
  call, so that the string and its length are the last two arguments in the
  corresponding C-function.

  For the definition of the open flag in the call "openf()", it is required
  to include the file "fcntl.h". On most operating systems, this file
  is in "/usr/include", but we have encountered systems where it is in
  "/usr/include/sys". To support both, you should compile this source-code
  with the additional search path specified: -I/usr/include/sys.

  The syntax of the Fortran calls are provided below.

-------
 2-Dec-1992 - fastio.c: John L. Lillibridge, NOAA/NOS/OES Geosciences Lab
11-Nov-1993 - Addition of seekf and warning by Remko Scharroo, DUT/SSR&T
17-May-1994 - Nice manual
 7-Jan-1996 - Underscores through UNDERSCOREAFTER. Warning also without
              carriage return.
12-Oct-2001 - Added line reading and writing. Add read/write flag to IOCONST.
20-Dec-2001 - Improved manual, added O_TRUNC to the read/write flag.
-------

*&OPENF -- Open file for FASTIO
*+
*     FUNCTION OPENF (OFLAG, MODE, FNAME)
*     INTEGER OPENF, OFLAG, MODE
*     CHARACTER*(*) FNAME
*
* Opens file FNAME with open flag OFLAG and permission mode MODE.
*
* OFLAG specifies opening for read or write. Most practical values are:
* O_RDONLY                     = o'0000' = Open READ ONLY, file must EXIST.
* O_WRONLY | O_CREAT | O_TRUNC = o'1101' = Open WRITE ONLY. If file exist
*        truncate it; if file does not exist, create it.
*
* The actual values may, however, differ per system. For the actual values
* on your machine consult the include file "/usr/include/fcntl.h", or use
* the subroutine IOCONST described below.
*
* MODE is an octal number describing the permission mode, as used with
* the "chmod" command. The permission code is only used when a file is
* created, and ignored in any other case. The permission code to be given
* is the logical AND of MODE and !UMASK. In most practical cases
* MODE = o'666' can be used.
*
* OPENF returns the file descriptor for use in subsequent calls to
* readf, writef, or closef. If OPENF is negative, an error occurred while
* opening the file. Use PERRORF to report the error.
*
* Usage for reading an EXISTING file.
*
*     INTEGER R_FLAG, W_FLAG, FD, OPENF, STRING
*
*     CALL IOCONST (R_FLAG, W_FLAG)
*     FD = OPENF (R_FLAG, O'666', 'input_file')
*     CALL READF (FD, 4, STRING)
*     CALL CLOSEF (FD)
*
* Usage for writing a file which may or may not exist.
*
*     FD = OPENF (W_FLAG, O'666', 'output_file')
*     CALL WRITEF (FD, 4, STRING)
*
* Note: Standard input and standard output are NOT opened through OPENF.
* Simply assign FD=0 for standard input and FD=1 for standard output.
*
* Input argument:
*  OFLAG : Open flag (see above or IOCONST).
*  MODE  : Permission mode. MODE = O'666' is recommended.
*  FNAME : File name of the input or output file. Maximum 1024 chars.
* Returned value:
*  OPENF : Contains the file descriptor on return. If OPENF < 0
*                   an error occurred while opening the file.
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
*  NBYTE  : Number of bytes to be read.
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
*  NBYTE  : Number of bytes to be written
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
* "COMMENT: error message", where COMMENT is hold in the parameter in the
* call to PERRORF.
*
* Input argument:
*  COMMENT : Additional comment to the error message
*-

*&IOCONST -- I/O flags for FASTIO routines
*+
*     SUBROUTINE IOCONST (R_FLAG, W_FLAG)
*     INTEGER R_FLAG, W_FLAG
*
* Returns I/O flags for read and write by using constants in
* "/usr/include/fcntl.h":
*
*      R_FLAG = O_RDONLY
*      W_FLAG = O_WRONLY + O_CREAT + O_EXCL
*     RW_FLAG = O_RDWR + O_CREAT + O_TRUNC
*
* Output arguments:
*  R_FLAG  : Equivalent to Fortran 'status=old'
*            Opens for reading. File must exist.
*  W_FLAG  : Equivalent to Fortran 'status=new'
*            Creates file for writing. File should not yet exist.
*  RW_FLAG : Equivalent to Fortran 'status=unknown'
*            Opens for writing. Creates file if it does not exist.
*            Truncates file to last written byte.
*-

*&WARNING -- Write string to standard error
*+
*     SUBROUTINE WARNING (STRING)
*     CHARACTER*(*) STRING
*
* Writes STRING to standard error using low-level I/O routines. Replaces
*     WRITE (0,'(a)') STRING  and
*     WRITE (2,'(a)') STRING.
* When STRING ends on /$ (slash dollar), the carriage return
* will not be parsed after the contents of STRING.
*
* Input argument:
*  STRING : Warning prompted to standard error
*-
*=
*/

#define stringlength 1024

#include <stdio.h>
#include <sys/types.h>
#include <fcntl.h>

#ifdef CAPITALS
#define openf OPENF
#define closef CLOSEF
#define readf READF
#define writef WRITEF
#define perrorf PERRORF
#define ioconst IOCONST
#define seekf SEEKF
#define warning WARNING
#endif

#ifdef UNDERSCOREBEFORE
#define openf _openf
#define closef _closef
#define readf _readf
#define writef _writef
#define perrorf _perrorf
#define ioconst _ioconst
#define seekf _seekf
#define warning _warning
#endif

#ifdef UNDERSCOREAFTER
#define openf openf_
#define closef closef_
#define readf readf_
#define writef writef_
#define perrorf perrorf_
#define ioconst ioconst_
#define seekf seekf_
#define warning warning_
#endif
/*****************************************************************************/

long openf(oflag,mode,fname,slen)
char *fname;
int slen, *oflag, *mode;
{
  char temp[stringlength];              /* Scratch string...could malloc() */
  int  len = slen - 1;                  /* String length variable */

  strncpy(temp,fname,slen);             /* Use copy of F77 string */

/* Strip off trailing blanks. */

  while (temp[len] == ' ' && len >= 0)
       len--;

  temp[len+1] = '\0';                  /* Null terminate string for C call */

  return(open(temp,*oflag,*mode));
}

/*****************************************************************************/

long closef(fd)
int *fd;
{
  if(*fd == 0){
    return(fclose(stdin));
  }else if(*fd == 1){
    return(fclose(stdout));
  }else{
    return(close(*fd));
  }
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
      fgets(buf,slen,stdin);	/* null-terminated, WITH <LF> char returned */
      return(strlen(buf));
    }else{			/* read nbytes of binary data */
      return(fread(buf,1,*nbytes,stdin));  /* Number of 1-byte items=*nbytes */
    }
  }else{			/* read from file descriptor *fd */
    if(*nbytes == 0){		/* read ascii "lines" */
      for(i=0; i<slen; i++){
        read(*fd,&temp[i],1);
        if(temp[i]==10)break;
      }
      i++;
      temp[i]='\0';		/* null-terminated, WITH <LF> char returned */
      strncpy(buf,temp,i);
      return(i);
     }else{			/* read nbytes of binary data */
      return(read(*fd,buf,*nbytes));
     }
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
      for(i=0; i<slen; i++){
        if(temp[i]==10)break;
      }
      i++; temp[i]='\0';
      return(fputs(temp,stdout));
    }else{			/* write nbytes of binary data */
      return(fwrite(buf,1,*nbytes,stdout)); /* Number of items = *nbytes */
    }
  }else{			/* write to file descriptor *fd */
    if(*nbytes == 0){		/* nbytes=0 signals write ascii "lines" */
      strncpy(temp,buf,slen);
      for(i=0; i<slen; i++){
        if(temp[i]==10)break;
      }
      i++; temp[i]='\0';
      return(write(*fd,temp,i));
    }else{			/* write nbytes of binary data */
      return(write(*fd,buf,*nbytes));
    }
  }
}

/*****************************************************************************/

void perrorf(string,slen)
char *string;
int slen;
{
  char temp[stringlength];              /* Scratch string...could malloc() */

  int len = slen - 1;

  strncpy(temp,string,slen);

/* Strip off trailing blanks. */

  while (temp[len] == ' ' && len >= 0)
       len--;

  temp[len+1] = '\0';                  /* Null terminate string for C call */

  perror(temp);
  return;
}

/*****************************************************************************/


void ioconst(r_flag,w_flag,rw_flag)
int *r_flag, *w_flag, *rw_flag;
{
  *r_flag = O_RDONLY;
  *w_flag = O_WRONLY | O_EXCL | O_CREAT;
  *rw_flag= O_RDWR;
  *rw_flag= O_RDWR | O_TRUNC | O_CREAT;
  return;
}

/*****************************************************************************/

long seekf(fd,nbytes,whence)
int *fd, *nbytes, *whence;
{
  if(*fd == 0){
    return(fseek(stdin,*nbytes,*whence));
  }else{
    return(lseek(*fd,*nbytes,*whence));
  }
}

/*****************************************************************************/

void warning(text,text_len)
char *text;
int text_len;
{
  int i;
  int last = text_len;

  while (last && text[last - 1] == ' ') 
       last--;
  
  if (last > 1 && text[last-2] == '/' && text[last-1] == '$')
  {
    for (i = 0; i < last-2; i++) { putc(text[i], stderr); }
  } else {
    for (i = 0; i < last; i++) { putc(text[i], stderr); }
    putc('\n', stderr);
  }
  fflush(stderr);
}
