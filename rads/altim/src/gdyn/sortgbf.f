      PROGRAM SORTGBF

C******************************************************************************
C  This program sorts the records of a data file in Geodyn Binary Format, in  *
C  acsending time order. The sorting is done by reading the data records and  *
C  puting them in the array IMASS, the times are written to the array TIMES,  *
C  which is then sorted. The indexes which connects the times with the data   *
C  records are kept in the array INDEX. The data records from the array IMASS *
C  are written to an output file using the INDEX array. When the number of    *
C  input data records is less or equal to the parameter MAXRES the output     *
C  file is the final output file. However, when there are more input data     *
C  records then MAXRES, the data records are first written to a scratch file, *
C  and the process is repeated until there are no more input data records     *
C  left or the maximum number of scratch files (MAXFIL) is used. After that   *
C  the scratch files are written to the final output file. The number of      *
C  input data records is limited to the product of the parameters MAXREC and  *
C  MAXFIL.                                                                    *
C                                                                             *
C  Input file:                                                                *
C     Unit 10 (gbf.in) : Data file in Geodyn Binary Format.                   *
C                                                                             *
C  Scratch files:                                                             *
C     Unit 30-.. : Scratch files (maximum number of files: MAXFIL).           *
C                                                                             *
C  Output file:                                                               *
C     Unit 20 (gbf.out) : Sorted data file in Geodyn Binary Format.           *
C                                                                             *
C  The input and output file names can be changed from the default by giving  *
C  the names as arguments: sortgbf input-filename -o output-filename.         *
C                                                                             *
C                                                                             *
C  Created by G.J.Mets, March 30, 1994                                        *
C                                                                             *
C******************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER ( MAXREC = 100000 , MAXFIL = 10 )

C  Parameter MAXREC controls the amount of memory needed for running the 
C  program. Memory needed is MAXREC * 80 bytes.
C  Parameter MAXFIL is the maximum number of scratch files. This number
C  multiplied with MAXREC gives the maximum number of oberservations
C  this program can sort.

      CHARACTER *80 ARG, FNIN, FNOUT
      REAL TARRAY(2), DUM
      LOGICAL LIN

      DIMENSION TMAX(MAXFIL), TMIN(MAXFIL), IDAT(17,MAXFIL),
     $          IFSCR(MAXFIL), TIME(MAXFIL)
      DIMENSION IDATA(17) , TIMES(MAXREC) , INDEX(MAXREC)
      DIMENSION IMASS(17,MAXREC)

      EQUIVALENCE ( IDATA(5) , MJD ) , ( IDATA(6) , TFRAC )

      DATA ARG, FNIN, FNOUT / ' ', 'gbf.in', 'gbf.out' /

C  Initialisation:
C     Unit numbers of input (10) and output (6,20) files.

      IF6    =  6
      IFIN   = 10
      IFOUT  = 20

C     Counters for: number of records NREC (large arrays)
C                   number of scratch files NSFIL
C                   number of observations NOBS

      NSFIL  = 0
      NREC   = 0
      NOBS   = 0

C  End initialisation.

C  Read the arguments if available.
      NRARG = IARGC()
      I = 0
 5    I = I + 1
      IF ( I .LE. NRARG ) THEN
        CALL GETARG ( I, ARG )
        IF ( ARG(1:2) .EQ. '-o' ) THEN
          I = I + 1
          CALL GETARG ( I, FNOUT )
        ELSE
          FNIN = ARG
        ENDIF
        GOTO 5
      ENDIF

C  Open input file: IFIN (default name is: gbf.in).

      INQUIRE ( FILE = FNIN, EXIST = LIN )
      IF ( .NOT. LIN ) THEN
        WRITE ( IF6, * ) 'ERROR: Input file is not present.'
        WRITE ( IF6, * ) '       Execution stops in SORTGBF.'
        STOP 99
      ENDIF
      OPEN ( IFIN, FILE = FNIN, FORM = 'UNFORMATTED' )

C  Reading observation-records from input file

 10   READ ( IFIN , END = 100 ) IDATA

      NREC = NREC + 1
      NOBS = NOBS + 1

C  If number of data records (NREC) is larger than MAXREC, sort the TIMES array
C  in ascending time order and write the sorted records to a scratch file.

      IF ( NREC .GT. MAXREC ) THEN
         NSFIL = NSFIL + 1

C  If number of scratch files is larger than MAXFIL, the number of observations
C  is too large and the program will stop.

         IF ( NSFIL .GT. MAXFIL ) CALL ERROR (1)
         IFSCR(NSFIL) = 29 + NSFIL

         OPEN ( IFSCR(NSFIL), STATUS = 'SCRATCH', FORM = 'UNFORMATTED' )
         CALL SORTIM ( MAXREC, TIMES, INDEX, IMASS, IFSCR(NSFIL),
     $                 TMAX(NSFIL), TMIN(NSFIL) )

C  Reset number of data records for the array IMASS

         NREC = 1

      ENDIF

C  Write the input data record to the array IMASS and write the
C  time and index to the arrays TIMES and INDEX.

      DO 50 J = 1 , 17
         IMASS(J,NREC) = IDATA(J)
 50   CONTINUE 
      TIMES(NREC) = DBLE(MJD) + TFRAC
      INDEX(NREC) = NREC

      GOTO 10

C  Reading of the observation records is completed. 

 100  CLOSE ( IFIN )

C  If the number of scratch files is zero, sort the TIMES array and
C  write the sorted records directly to the output file IFOUT. If not zero,
C  sort the data records and write them to a new scratch file.

      IF ( NSFIL .EQ. 0 ) THEN

         OPEN ( IFOUT, FILE = FNOUT, FORM = 'UNFORMATTED' )
         CALL SORTIM ( NREC, TIMES, INDEX, IMASS, IFOUT, DUM1, DUM2 )
         GOTO 500

      ELSE

         NSFIL = NSFIL + 1
         IF ( NSFIL .GT. MAXFIL ) CALL ERROR (1)
         IFSCR(NSFIL) = 29 + NSFIL

         OPEN ( IFSCR(NSFIL), STATUS = 'SCRATCH', FORM = 'UNFORMATTED' )
         CALL SORTIM ( NREC, TIMES, INDEX, IMASS, IFSCR(NSFIL),
     $                 TMAX(NSFIL), TMIN(NSFIL) )

      ENDIF

C  If scratch files were used, write the sorted observation records to the
C  output file (IFOUT) by taking each time the earliest observation time.

C  First sort the scratch files in ascending order looking at the latest
C  observation time of each file.

      CALL SORTFL ( IFSCR, NSFIL, TMAX, TMIN )

C  Read the first record of all scratch files and determine the times of the 
C  observation records.

      DO 120 I = 1, NSFIL
         REWIND ( IFSCR(I) )
         READ ( IFSCR(I) ) IDATA
         TIME(I) = DBLE ( MJD ) + TFRAC
         DO 110 J = 1 , 17
            IDAT(J,I) = IDATA(J)
 110     CONTINUE
 120  CONTINUE 

      OPEN ( IFOUT, FILE = FNOUT, FORM = 'UNFORMATTED' )

C  Determine which data record comes first and write this record to the output
C  file (IFOUT). Read new record from the scratch file in question.

      ISTART = 1
 130  ISTART = ISTART + 1
 140  IMIN   = ISTART - 1

      DO 150 I = ISTART , NSFIL
         IF ( TIME(IMIN) .GT. TIME(I) ) IMIN = I
 150  CONTINUE 

      WRITE ( IFOUT ) ( IDAT(J,IMIN) , J = 1 , 17 )
      READ ( IFSCR(IMIN) , END = 160 ) IDATA
      TIME(IMIN) = DBLE ( MJD ) + TFRAC
      DO 155 J = 1 , 17
         IDAT(J,IMIN) = IDATA(J)
 155  CONTINUE
      GOTO 140

C  If all records from a scratch file are written, close the file.

 160  CLOSE ( IFSCR(IMIN) )
      IF ( IMIN .EQ. NSFIL ) GOTO 500
      GOTO 130

 500  CONTINUE 
      CLOSE ( IFOUT )

      WRITE ( IF6 , 800 ) NOBS
 800  FORMAT ( I7, ' records are sorted.' )
      WRITE ( IF6 , * ) 'Normal end of program SORTGBF.'

      END
*******************************************************************************
      SUBROUTINE SORTIM ( NREC, TIMES, INDEX, IMASS, IFOUT, TMAX, TMIN )

C  This subroutine sorts the array TIMES in ascending order. The array INDEX
C  is used to keep the connection between the times in TIMES and the data
C  records in the array IMASS. After sorting the array IMASS is written to
C  IFOUT.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TIMES(NREC), INDEX(NREC), IDATA(17), IMASS(17,NREC)

      DO 100 I = 2 , NREC
         IND = I
         IF ( TIMES(IND) .LT. TIMES(IND-1) ) THEN
            THULP = TIMES(IND)
            IHULP = INDEX(IND)
 10         TIMES(IND) = TIMES(IND-1)
            INDEX(IND) = INDEX(IND-1)
            IND = IND - 1
            IF ( IND .EQ. 1 ) GOTO 20
            IF ( THULP .LT. TIMES(IND-1) ) GOTO 10
 20         TIMES(IND) = THULP
            INDEX(IND) = IHULP
         ENDIF
 100  CONTINUE 

C  Writing IMASS to IFOUT

      DO 110 I = 1 , NREC
         DO 105 J = 1 , 17
            IDATA(J) = IMASS(J,INDEX(I))
 105     CONTINUE 
         WRITE ( IFOUT ) IDATA
 110  CONTINUE

      TMAX = TIMES(NREC)
      TMIN = TIMES(1)

      END
*******************************************************************************
      SUBROUTINE SORTFL ( IFSCR , NSFIL , TMAX , TMIN )

C  This subroutine sorts the sequence of the scratch files so that the times 
C  in the array TMAX are in ascending order.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IFSCR(NSFIL) , TMAX(NSFIL) , TMIN(NSFIL)

      DO 100 I = 2 , NSFIL
         IND = I
         IF ( TMAX(IND) .LT. TMAX(IND-1) ) THEN
            THULP1 = TMAX(IND)
            THULP2 = TMIN(IND)
            IHULP  = IFSCR(IND)
 10         TMAX(IND) = TMAX(IND-1)
            TMIN(IND) = TMIN(IND-1)
            IFSCR(IND) = IFSCR(IND-1)
            IND = IND - 1
            IF ( IND .EQ. 1 ) GOTO 20
            IF ( THULP1 .LT. TMAX(IND-1) ) GOTO 10
 20         TMAX(IND) = THULP1
            TMIN(IND) = THULP2
            IFSCR(IND) = IHULP
         ENDIF
 100  CONTINUE

      END
*******************************************************************************
      SUBROUTINE ERROR ( I )

C This subroutine prints error messages and stops the execution of the program.

      IF6  = 6

      GOTO ( 10 ) I

 10   WRITE ( IF6 , 15 )
 15   FORMAT ( ' ERROR: Total number of observations is to large.'/
     $         8X,'Increase the value of the parameters MAXREC/MAXFIL.'/
     $         8X,'Execution stops in programs sort/gbf.' )
      STOP 99

      END
