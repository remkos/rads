      PROGRAM SUMMG01
************************************************************************
*
*     This program gives a summary of a GEODYN run. The program has
*     an option for the printout of:
*     - the inner iteration(s) in the final outer iteration;
*     - the final station coordinates (earth fixed rectangular and
*       geodetic);
*     - the final solutions for the earth rotation parameters;
*     - the back-substitution deltas of the arc parameters after the
*       final outer iteration.
*
*     Irrespective of the input parameter, always printed are the
*     current changes in:
*     - the rectangular station coordinates:
*     - the earth rotation parameter solutions.
*
*     The program can handle multiple data arcs.
*
*     Input:
*     FT10: a file with the Geodyn printout (preferably GEODYN I 8210)
*     FT11: the corresponding residual file.
*     SYSIN: a parameter which determines whether (0) or not (1)
*            parts of the Geodyn printout are to be printed.
*
*     Output:
*     none.
*
*     N.B.
*     For running this program at the CONVEX computer remove the
*     C's in front of the OPEN ( ...., FILE =....) statements.
*
*
*     R. Noomen,   November 7, 1989.
*     R. Noomen,   January 15, 1990.
*     R. Noomen,    August 20, 1990.
*     W.J. van Gaalen, June 3, 1991.
*                      ( Modified for processing GEODYN II printout )
*     G.J. Mets,    January 3, 1992.
*             (The program can now also handle GEODYN II printout with
*              only inner iterations. From now on the drag tables are
*              skipped in the printout of the summary)
*
*     W.J. van Gaalen   October 14, 1992.
*             (Modified for processing of multi arc - multi satellite
*             of GEODYN 9204 printout.)
*
*     R. Noomen         November 1, 1992.
*             (Modified for processing of multi arc - multi satellite
*             of GEODYN 9208 printout.)
*
*     W.J. van Gaalen   November 11, 1992.
*             (Modified for multiple page printing of the residual
*             summary by station.)
*
*
**********************************************************************
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
      CHARACTER * 133 C133
      CHARACTER * 10 ARG
      LOGICAL SHORT
C
      PARAMETER ( NDIM1 = 100 , NDIM2 = 10 )
      DIMENSION NITER(NDIM2),NPAGE(NDIM2)
      DIMENSION STDERP(NDIM1,NDIM2,3) , XERP(NDIM1,NDIM2,3) ,
     .          XAPERP(NDIM1,3)
      DIMENSION T1ERP(NDIM1) , T2ERP(NDIM1)
      DIMENSION STDSTA(NDIM1,NDIM2,3) , XSTA(NDIM1,NDIM2,3) ,
     .          XAPSTA(NDIM1,3)
      DIMENSION ISTAT(NDIM1)
      DIMENSION DELTA(3) , SIGMA(3)
C
      DATA IFIN , IFOUT , IF5 , IF6 / 10 , 6 , 5 , 6 /
      DATA NARCS , NITER1 , IGDYN , NSTAT , NERP , JPARM , IPRINT
     .  /  6 * 0 , 2 /
      DATA NARCHL , NITHLP , IOUTHL
     .  /  3 * 0  /
C
      DPI = 4.0D0 * ATAN ( 1.0D0 )
      SHORT = .FALSE.
C
C  For Convex: remove C's of Comment in next statements
      IFIN = 5
      IPRINT=0
      narg=iargc()
      do 5 i=1,narg
	 call getarg(i,arg)
         IF (ARG.EQ.'-l') IPRINT=0
         IF (ARG.EQ.'-s') IPRINT=1
	 if (arg.eq.'-h') then
	    write (if6, *) 'usage: summg01 [ -l | -s ]'
	    goto 9999
	 endif
    5 continue

C     OPEN ( IFIN  , FILE  = 'printout' )
C     OPEN ( IFOUT , FILE  =  'summary' )
C     OPEN ( IF5   , FILE  =    'sysin' )
C     OPEN ( IF6   , FILE  =   'output' )
C
C     READ   ( IF5 , * , END = 998 ) IPRINT
      IF ( IPRINT .NE. 0 .AND. IPRINT .NE. 1 ) THEN
        WRITE  ( IF6 , * ) 'incorrect value for input parameter ' ,
     .    'IPRINT:' , IPRINT , '. stop.'
        GOTO 9999
      END IF
C
      REWIND IFIN
C
C  determine the GEODYN version that is used,
C  the number of data arcs, the number of iterations
C  and the ids of the stations that are to be estimated.
      REWIND IFIN
  10  READ   ( IFIN , 20 , END = 600 ) C133
  20  FORMAT ( 1A133 )
      IF ( C133(1:4) .EQ. '1X  ' ) THEN
        CALL READI4 ( IF6 , C133 , 41 , 44 , IGDYN )
        IF ( IGDYN .NE. 8202 .AND. IGDYN .NE. 8210 ) THEN
          WRITE  ( IF6 , * ) 'Unknown GEODYN version' , IGDYN ,
     .      '. Stop.'
          GOTO 9999
        END IF
        IF ( IGDYN .EQ. 8202 ) WRITE  ( IF6 , * )
     .    'WARNING: the proper working of this program is not tested ',
     .    'thoroughly for GEODYN I version 8202. Execution continues.'
        IF ( IGDYN .EQ. 8202 ) THEN
          IPLUS = 0
        ELSE
          IPLUS = 1
        END IF
      END IF
      IF ( C133(2:11) .EQ. 'GEODYN IIE' ) THEN
        CALL READI4 ( IF6 , C133 , 23 , 26 , IGDYN )
        IF ( IGDYN .NE. 8901 .AND.
     .       IGDYN .NE. 9204 .AND.
     .       IGDYN .NE. 9208 ) THEN
          WRITE  ( IF6 , * ) 'Unknown GEODYN IIE version' , IGDYN ,
     .      '. Stop.'
          GOTO 9999
        END IF
        IPLUS = 1
      END IF
C
      IF ( IGDYN .EQ. 8901 .OR.
     .     IGDYN .EQ. 9204 .OR.
     .     IGDYN .EQ. 9208 ) THEN
        IF (C133(13:28) .EQ. 'RESIDUAL SUMMARY' ) THEN
   25     READ (C133,30) NARCS,NIT,IOUTER
          IF (NARCS.EQ.NARCHL.AND.NIT.EQ.NITHLP.AND.IOUTER.EQ.IOUTHL)
     .      THEN
            NPAGE(NARCS) = NPAGE(NARCS) + 1
          ELSE
            NPAGE(NARCS) = 1
          END IF
          NARCHL = NARCS
          NITHLP = NIT
          IOUTHL = IOUTER
          NITER(NARCS) = NIT
   30     FORMAT (56X,1I3,16X,1I3,20X,1I2)
          IF ( NARCS .GT. NDIM2 ) THEN
            WRITE  ( IF6 , * )
     .      'The number of arcs exceeds',NDIM2 ,'. STOP.'
            GOTO 9999
          END IF
          IF ( IOUTER .GT. NDIM2 ) THEN
            WRITE  ( IF6 , * )
     .      'The number of outer iterations exceeds' , NDIM2 , '. STOP.'
            GOTO 9999
          END IF
          GOTO 400
        END IF
      ELSE
        IF ( C133(2:7) .EQ. 'DATA  ' ) NARCS = NARCS + 1
        IF ( C133(2:7) .EQ. 'ENDSTA' ) NARCS = NARCS + 1
        IF ( C133(2:7) .EQ. 'ENDGLB' ) NARCS = NARCS + 1
        IF ( C133(1:4) .EQ. '1V  '   ) NITER1 = NITER1 + 1
        IF ( C133(37:52) .EQ. 'STATION POSITION'
     .    .OR. C133(67:81) .EQ. 'POLE AND A1-UT1'
     .    .OR. C133(64:78) .EQ. 'POLE AND A1-UT1' ) THEN
          IOUTER = NITER1 / NARCS
          IF ( IOUTER .GT. NDIM2 ) THEN
            WRITE  ( IF6 , * )
     .      'THE NUMBER OF OUTER ITERATIONS EXCEEDS' , NDIM2 , '. STOP.'
            GOTO 9999
          END IF
          GOTO 300
        END IF
      END IF
      GOTO 10
C----------------------------------------------------------------------
C====================    GEODYN I    ==================================
C----------------------------------------------------------------------
C
 300  IF ( C133(37:52) .EQ. 'STATION POSITION' ) THEN
 310    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(6:13)  .EQ. 'GEODETIC' ) GOTO 10
        IF ( C133(18:25) .EQ. 'CYLINDRI' ) GOTO 10
C
        IF ( C133(2:9) .EQ. 'A PRIORI' .AND. IOUTER .EQ. 1 ) THEN
          IF ( IGDYN .EQ. 8202 ) THEN
            READ   ( C133 , 320 ) ITEST
 320        FORMAT ( 21X , 1I4 )
          ELSE
            READ   ( C133 , 330 ) ITEST
 330        FORMAT ( 22X , 1I4 )
          END IF
          NSTAT = NSTAT + 1
          IF ( NSTAT .GT. NDIM1 ) THEN
            WRITE  ( IF6 , * ) 'The number of estimated stations' ,
     .        'exceeds' , NDIM1 , '. Stop.'
            GOTO 9999
          END IF
          ISTAT(NSTAT) = ITEST
          CALL READR8
     .      ( IF6 , C133 , 27 + IPLUS , 40 + IPLUS , XAPSTA(NSTAT,1) )
          CALL READR8
     .      ( IF6 , C133 , 41 + IPLUS , 54 + IPLUS , XAPSTA(NSTAT,2) )
          CALL READR8
     .      ( IF6 , C133 , 55 + IPLUS , 68 + IPLUS , XAPSTA(NSTAT,3) )
        END IF
C
        IF ( C133(2:9) .EQ. 'ADJUSTED' ) THEN
          IF ( IGDYN .EQ. 8202 ) THEN
            READ   ( C133 , 320 ) ITEST
          ELSE
            READ   ( C133 , 330 ) ITEST
          END IF
          DO 340 IJ = 1 , NSTAT
            IF ( ITEST .EQ. ISTAT(IJ) ) THEN
              CALL READR8 ( IF6 , C133 , 27 + IPLUS , 40 + IPLUS ,
     .          XSTA(IJ,IOUTER,1) )
              CALL READR8 ( IF6 , C133 , 41 + IPLUS , 54 + IPLUS ,
     .          XSTA(IJ,IOUTER,2) )
              CALL READR8 ( IF6 , C133 , 55 + IPLUS , 68 + IPLUS ,
     .          XSTA(IJ,IOUTER,3) )
              CALL READR8 ( IF6 , C133 , 71 + IPLUS , 80 + IPLUS ,
     .          STDSTA(IJ,IOUTER,1) )
              CALL READR8 ( IF6 , C133 , 81 + IPLUS , 90 + IPLUS ,
     .          STDSTA(IJ,IOUTER,2) )
              CALL READR8 ( IF6 , C133 , 91 + IPLUS , 100 + IPLUS ,
     .          STDSTA(IJ,IOUTER,3) )
              GOTO 310
            END IF
 340      CONTINUE
        END IF
        GOTO 310
      END IF
C
C  readout of a priori values for earth rotation parameters
C
      IF ( C133(64:78) .EQ. 'POLE AND A1-UT1' ) THEN
 350    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(58:76) .EQ. 'TRACKING COMPLEMENT' ) GOTO  10
        IF ( C133(58:76) .EQ. '                   ' ) GOTO 350
        IF ( C133(65:70) .EQ. 'A1-UT1' ) GOTO 350
        IF ( C133(65:70) .EQ. ' (SEC)' ) GOTO 350
        NERP = NERP + 1
        IF ( NERP .GT. NDIM1 ) THEN
          WRITE  ( IF6 , * )
     .      'The number of estimated ERPs exceeds' , NDIM1 , '. Stop.'
          GOTO 9999
        END IF
        CALL READR8 ( IF6 , C133 ,  3 , 12 , T1ERP(NERP) )
        CALL READR8 ( IF6 , C133 , 18 , 27 , T2ERP(NERP) )
        CALL READR8 ( IF6 , C133 , 32 , 42 , XAPERP(NERP,1) )
        CALL READR8 ( IF6 , C133 , 47 , 57 , XAPERP(NERP,2) )
        CALL READR8 ( IF6 , C133 , 64 , 70 , XAPERP(NERP,3) )
        GOTO 350
      END IF
C
C  readout of earth rotation parameter solutions
C
      IF ( C133(67:81) .EQ. 'POLE AND A1-UT1' ) THEN
        DO 370 IERP = 1 , NERP
 360      READ   ( IFIN , 20 , END = 999 ) C133
          IF ( C133(58:76) .EQ. '                   ' ) GOTO 360
          IF ( C133(65:70) .EQ. 'A1-UT1' ) GOTO 360
          IF ( C133(65:70) .EQ. ' (SEC)' ) GOTO 360
          CALL READR8 ( IF6 , C133 ,  3 , 12 , TEST1 )
          CALL READR8 ( IF6 , C133 , 18 , 27 , TEST2 )
          TTEST = ABS ( TEST1 - T1ERP(IERP) )
     .      + ABS ( TEST2 - T2ERP(IERP) )
          IF ( TTEST .GT. 1.0D-5 ) THEN
            WRITE  ( IF6 , * ) 'Incorrect order of ERPs. Stop.'
            GOTO 9999
          END IF
          CALL READR8 ( IF6 , C133 ,  32 ,  42 ,   XERP(IERP,IOUTER,1) )
          CALL READR8 ( IF6 , C133 ,  47 ,  57 ,   XERP(IERP,IOUTER,2) )
          CALL READR8 ( IF6 , C133 ,  64 ,  70 ,   XERP(IERP,IOUTER,3) )
          CALL READR8 ( IF6 , C133 ,  89 ,  98 , STDERP(IERP,IOUTER,1) )
          CALL READR8 ( IF6 , C133 , 104 , 113 , STDERP(IERP,IOUTER,2) )
          CALL READR8 ( IF6 , C133 , 119 , 128 , STDERP(IERP,IOUTER,3) )
 370    CONTINUE
      END IF
C
      GOTO 10
C----------------------------------------------------------------------
C====================    GEODYN II   ==================================
C----------------------------------------------------------------------
C
 400  IF ( C133(12:28) .EQ. 'FIXED RECTANGULAR' ) THEN
 410    READ   ( IFIN , 20 , END = 600 ) C133
        IF ( C133(6:13) .EQ. 'GEODETIC' ) GOTO 10
        IF ( C133(18:28) .EQ. 'CYLINDRICAL' ) GOTO 10
C
        IF ( C133(2:9) .EQ. ' APRIORI' .AND. IOUTER .EQ. 1 ) THEN
          READ   ( C133 , 330 ) ITEST
          NSTAT = NSTAT + 1
          IF ( NSTAT .GT. NDIM1 ) THEN
            WRITE  ( IF6 , * ) 'The number of estimated stations' ,
     .        'exceeds' , NDIM1 , '. Stop.'
            GOTO 9999
          END IF
          ISTAT(NSTAT) = ITEST
          CALL READR8
     .      ( IF6 , C133 , 27 + IPLUS , 40 + IPLUS , XAPSTA(NSTAT,1) )
          CALL READR8
     .      ( IF6 , C133 , 41 + IPLUS , 54 + IPLUS , XAPSTA(NSTAT,2) )
          CALL READR8
     .      ( IF6 , C133 , 55 + IPLUS , 68 + IPLUS , XAPSTA(NSTAT,3) )
        END IF
C
        IF ( C133(2:9) .EQ. 'ADJUSTED' ) THEN
            READ   ( C133 , 330 ) ITEST
          DO 440 IJ = 1 , NSTAT
            IF ( ITEST .EQ. ISTAT(IJ) ) THEN
              CALL READR8 ( IF6 , C133 , 27 + IPLUS , 40 + IPLUS ,
     .          XSTA(IJ,IOUTER,1) )
              CALL READR8 ( IF6 , C133 , 41 + IPLUS , 54 + IPLUS ,
     .          XSTA(IJ,IOUTER,2) )
              CALL READR8 ( IF6 , C133 , 55 + IPLUS , 68 + IPLUS ,
     .          XSTA(IJ,IOUTER,3) )
              CALL READR8 ( IF6 , C133 , 71 + IPLUS , 80 + IPLUS ,
     .          STDSTA(IJ,IOUTER,1) )
              CALL READR8 ( IF6 , C133 , 81 + IPLUS , 90 + IPLUS ,
     .          STDSTA(IJ,IOUTER,2) )
              CALL READR8 ( IF6 , C133 , 91 + IPLUS , 100 + IPLUS ,
     .          STDSTA(IJ,IOUTER,3) )
              GOTO 410
            END IF
 440      CONTINUE
        END IF
        GOTO 410
      END IF
C
  500 READ (IFIN,20,END=600) C133
      IF (C133(2:29).EQ.'PARAMETER ADJUSTMENT SUMMARY') GOTO 500
      NERP=0
      JPARM=0
C
C  Read 'POLE X' and Time of Pole solutions
  510 READ (IFIN,20,END=600) C133
      IF (C133(1:5).EQ.'0POLE'.OR.C133(2:5).EQ.'A1 U'
     .  .OR.C133(2:5).EQ.'STA ') JPARM=JPARM+1
      IF (JPARM.GT.500) THEN
        WRITE (IF6,*) 'THE NUMBER OF ESTIMATED COMMON PARAMETERS ',
     .    'EXCEEDS 500. STOP.'
        GOTO 9999
      END IF
      IF (C133(1:7).EQ.'0POLE X') THEN
        NERP = NERP + 1
        IF ( NERP .GT. NDIM1 ) THEN
          WRITE  ( IF6 , * )
     .      'The number of estimated ERPs exceeds' , NDIM1 , '. Stop.'
          GOTO 9999
        END IF
        CALL READR8(IF6,C133,19,41,XAPERP(NERP,1))
        CALL READI4(IF6,C133,12,17,ITEST)
        MJD=MDATE(2,ITEST)
        T1ERP(NERP)=DBLE(MJD)-2.5D0
        T2ERP(NERP)=DBLE(MJD)+2.5D0
        MJD1=INT(T1ERP(NERP)+1.0D-10)
        MJD2=INT(T2ERP(NERP)+1.0D-10)
        DT1=T1ERP(NERP)-DBLE(MJD1)
        DT2=T2ERP(NERP)-DBLE(MJD2)
        IYMD1=MDATE(1,MJD1)
        IYMD2=MDATE(1,MJD2)
        T1ERP(NERP)=DBLE(IYMD1)+DT1
        T2ERP(NERP)=DBLE(IYMD2)+DT2
        DO 530 J=1,2
  520   READ (IFIN,20,END=600) C133
        IF (C133(112:118).EQ.'PAGE NO') GOTO 520
  530   CONTINUE
        CALL READR8(IF6,C133,19,41,XERP(NERP,IOUTER,1))
        CALL READR8(IF6,C133,63,77,STDERP(NERP,IOUTER,1))
C convert ERP(X) to arcsec
        XAPERP(NERP,1)=XAPERP(NERP,1)*3600.0D0*180.0D0/DPI
        XERP(NERP,IOUTER,1)=XERP(NERP,IOUTER,1)*3600.0D0*180.0D0/DPI
        STDERP(NERP,IOUTER,1)=STDERP(NERP,IOUTER,1)*3600.0D0*180.0D0/DPI
      END IF
C
C  Read 'POLE Y'
      IF (C133(1:7).EQ.'0POLE Y') THEN
        CALL READR8(IF6,C133,19,41,XAPERP(NERP,2))
        DO 545 J=1,2
  540   READ (IFIN,20,END=600) C133
        IF (C133(112:118).EQ.'PAGE NO') GOTO 540
  545   CONTINUE
        CALL READR8(IF6,C133,19,41,XERP(NERP,IOUTER,2))
        CALL READR8(IF6,C133,63,77,STDERP(NERP,IOUTER,2))
C convert ERP(Y) to arcsec
        XAPERP(NERP,2)=XAPERP(NERP,2)*3600.0D0*180.0D0/DPI
        XERP(NERP,IOUTER,2)=XERP(NERP,IOUTER,2)*3600.0D0*180.0D0/DPI
        STDERP(NERP,IOUTER,2)=STDERP(NERP,IOUTER,2)*3600.0D0*180.0D0/DPI
      END IF
C
C  Read 'A1 UT '
      IF (C133(2:7).EQ.'A1 UT ') THEN
        CALL READR8(IF6,C133,19,41,XAPERP(NERP,3) )
        DO 570 J=1,2
  560   READ (IFIN,20,END=600) C133
        IF (C133(112:118).EQ.'PAGE NO') GOTO 560
  570   CONTINUE
        CALL READR8(IF6,C133,19,41,XERP(NERP,IOUTER,3))
        CALL READR8(IF6,C133,63,77,STDERP(NERP,IOUTER,3))
      END IF
C
C     IF ( C133(12:28) .EQ. 'FIXED CYLINDRICAL' ) GOTO 10
      IF (C133(2:4).EQ.'C  ') GOTO 410
      IF (C133(13:28) .EQ. 'RESIDUAL SUMMARY' ) GOTO 25
C     IF (C133(2:4).EQ.'ARC') GOTO 10
      GOTO 510
C-----------------------------------------------------------------------
C
C  inventory ended
C
C  print interesting parts of GEODYN printout
 600  REWIND IFIN
C
      IF ( IGDYN .NE. 8901 .AND.
     .     IGDYN .NE. 9204 .AND.
     .     IGDYN .NE. 9208 ) THEN
C-----------------------------------------------------------------------
C  input deck Geodyn I
C-----------------------------------------------------------------------
 610    READ   ( IFIN , 20 , END = 999 ) C133
        IF (C133(2:5).NE.'X   ' ) GOTO 610
        call schrijf(iprint,ifout,c133)
 620    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(37:67) .NE. 'MULTI-ARC GEODYN RUN USING DATA' ) THEN
          call schrijf(iprint,ifout,c133)
          GOTO 620
        END IF
C       IF ( IPRINT .EQ. 0 ) THEN
C         WRITE  ( IFOUT , 630 )
C630      FORMAT ( '1' )
C       END IF
C
C  inner iteration(s) of final outer iteration
        NV = NITER1 - NARCS
        DO 640 I = 1 , NV
          CALL SURGE ( IFIN , 'V   ' , 'ZZZZ' , *999 , *999 )
 640    CONTINUE
 650    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(9:14) .NE. 'NSTEPS' ) GOTO 650
 660    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(3:11) .EQ. 'UNITS FOR' ) GOTO 660
        IF ( C133(2:10) .EQ. 'TRACEBACK' ) THEN
 665      READ ( IFIN , 20 , END = 999 ) C133
          IF (C133(2:20) .NE. 'STANDARD CORRECTIVE' ) GOTO 665
        END IF
        IF (C133(2:5).NE.'V   '.AND.C133(2:5).NE.'9   ') THEN
          call schrijf(iprint,ifout,c133)
          GOTO 660
        END IF
C
C  station coordinates in final outer iteration
 670    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(2:5) .NE. 'C   ' ) GOTO 670
        call schrijf(iprint,ifout,c133)
 680    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(18:28) .NE. 'CYLINDRICAL' ) THEN
          call schrijf(iprint,ifout,c133)
          GOTO 680
        END IF
C
C  polar motion and earth rotation in final outer iteration
        IF ( NERP .NE. 0 ) THEN
 690      READ   ( IFIN , 20 , END = 999 ) C133
          IF ( C133(2:5) .NE. 'Z   ' ) GOTO 690
          call schrijf(iprint,ifout,c133)
        END IF
 700    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(42:57) .NE. 'UPDATED ELEMENTS' ) THEN
          call schrijf(iprint,ifout,c133)
          GOTO 700
        END IF
C
C  back-substitution for all data arcs
        DO 730 I = 1 , NARCS
 710      READ   ( IFIN , 20 , END = 999 ) C133
          IF ( C133(46:54) .NE. 'BACK-SUBS' ) GOTO 710
          DO 715 J= 1,3
            call schrijf(iprint,ifout,c133)
            READ   ( IFIN , 20 , END = 999 ) C133
 715      CONTINUE
 720        READ   ( IFIN , 20 , END = 999 ) C133
          IF ( C133(8:24) .NE. 'CROSS CORRELATION' ) THEN
            IF ( IPRINT .EQ. 0 ) WRITE  (IFOUT,725) C133(8:50),
     .           C133(51:62),C133(63:80)
 725        FORMAT (1X,A43,17X,A12,15X,A18)
            GOTO 720
          END IF
 730    CONTINUE
C
      ELSE
C-----------------------------------------------------------------------
C  input deck  Geodyn II
C-----------------------------------------------------------------------
 740    READ   ( IFIN , 20 , END = 999 ) C133
        IF (C133(9:19).NE.'IIE VERSION' ) GOTO 740
        call schrijf(iprint,ifout,c133)
 750    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(2:7) .EQ. 'ETIDES' ) THEN
          IETID = IETID + 1
          GOTO 750
        END IF
        IF ( C133(2:7) .EQ. 'OTIDES' ) THEN
          IOTID = IOTID + 1
          GOTO 750
        END IF
        IF ( C133(2:10) .EQ. 'GEODYN II' ) GOTO 750
        IF ( C133 .EQ. ' ' ) GOTO 750
        IF (IETID.NE.0.AND.IOTID.NE.0.AND.IPRINT.EQ.0) THEN
          WRITE(IFOUT , 752 ) IETID
          WRITE(IFOUT , 753 ) IOTID
 752      FORMAT(I5,2X,'ETIDE CARDS SUPPRESSED')
 753      FORMAT(I5,2X,'OTIDE CARDS SUPPRESSED')
          IETID = 0
          IOTID = 0
        END IF
C       IF ( C133(58:65) .NE. 'TRACKING' ) THEN
C Modification by G.J.Mets January 3, 1992 (skipping drag tables)
        IF ( C133(58:65) .NE. 'TRACKING' .AND. C133(32:41) .NE.
     .       'NORMALIZED' ) THEN
C end of modification
          call schrijf(iprint,ifout,c133)
          GOTO 750
        END IF
C
C  Determine number of residual summaries to skip before final outer
C  iteration is reached
        NN = 0
        DO 770 I = 1 , NARCS
          NN = NITER(I)*NPAGE(I) + NN
 770    CONTINUE
        NN = NN*(IOUTER-1)
C
C  Skip all iterations in all outer iterations but the latest.
        DO 780 I = 1 , NN
 775      READ ( IFIN , 20 , END = 999 ) C133
          IF( C133(13:28) .NE. 'RESIDUAL SUMMARY' ) GOTO 775
 780    CONTINUE
C
        DO 799 I = 1 , NARCS
          IF (IOUTER.EQ.1.AND.JPARM.EQ.0.AND.NERP.EQ.0) SHORT=.TRUE.
C
C  Skip all inner iterations residual summary but the latest per arc.
          DO 784 J = 1 , (NITER(I)-1)*NPAGE(I)
 783        READ ( IFIN , 20 , END = 999 ) C133
            IF( C133(13:28) .NE. 'RESIDUAL SUMMARY' ) GOTO 783
 784      CONTINUE
C
C  Inner iteration(s) of final outer iteration.
 785      READ ( IFIN , 20 , END = 999 ) C133
          IF( C133(13:28) .NE. 'RESIDUAL SUMMARY' ) GOTO 785
          call schrijf(iprint,ifout,c133)
C
 790      READ   ( IFIN , 20 , END = 999 ) C133
          IF ( C133(51:62) .NE. 'OUTPUT UNITS' ) THEN
            call schrijf(iprint,ifout,c133)
            GOTO 790
          END IF
C
 796      READ   ( IFIN , 20 , END = 999 ) C133
          IF (SHORT) THEN
            IF ( C133(9:17) .NE. 'PARAMETER' ) GOTO 796
            call schrijf(iprint,ifout,c133)
            SHORT = .FALSE.
          ELSE
C           IF ( C133(42:60).EQ.'KEPLERIAN SATELLITE' ) GOTO 797
            IF ( C133(40:57).EQ.'SINGULAR KEPLERIAN'  ) GOTO 797
            IF ( C133(42:60).EQ.'CARTESIAN SATELLITE' ) GOTO 797
            IF(C133(2:5).EQ.'ARC '.AND.C133(9:15).EQ.'SUMMARY') GOTO 799
            GOTO 796
 797        call schrijf(iprint,ifout,c133)
          END IF
 798      READ   ( IFIN , 20 , END = 999 ) C133
C         IF (C133(2:5).EQ.'ARC '.AND.C133(9:15).EQ.'SUMMARY') GOTO 799
          IF (C133(46:49).NE.'PERG'.AND.C133(46:49).NE.'YDOT'.AND.
     .      C133(46:49).NE.'SICN') THEN
            call schrijf(iprint,ifout,c133)
            GOTO 798
          ELSE
            call schrijf(iprint,ifout,c133)
	    write (6,550)
            GOTO 796
          END IF
 799    CONTINUE
C
        IF (IOUTER.EQ.1.AND.JPARM.EQ.0.AND.NERP.EQ.0) GOTO 9999
C
C  station coordinates in final outer iteration
 800    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(2:5) .NE. 'C   ' ) GOTO 800
        call schrijf(iprint,ifout,c133)
 810    READ   ( IFIN , 20 , END = 999 ) C133
        IF ( C133(18:28) .NE. 'CYLINDRICAL' ) THEN
          call schrijf(iprint,ifout,c133)
          GOTO 810
        END IF
C
C  polar motion and earth rotation in final outer iteration
        IF ( NERP .NE. 0 .AND. IPRINT .EQ. 0) THEN
          WRITE (IFOUT,820)
     .     '1','ADJUSTED COORDINATES OF THE POLE AND A1-UT1 DIFFERENCE'
 820      FORMAT ( 1A1,36X,1A54)
          WRITE (IFOUT,822) 'BEGIN','END','X','Y','A1-UT1',
     .     'STANDARD DEVIATIONS'
 822      FORMAT (5X,A5,11X,A3,13X,A1,14X,A1,12X,A6,16X,A18)
          WRITE (IFOUT,824) '(SEC)','X','Y','A1-UT1'
 824      FORMAT (2(2X,'YYMMDD.DDD',3X),2(2X,' (ARC SEC)',3X),6X,A5,
     .      11X,A1,14X,A1,11X,A6)
          DO 827 I=1,NERP
           WRITE (IFOUT,825) T1ERP(I),T2ERP(I),(XERP(I,IOUTER,K),K=1,3),
     .          (STDERP(I,IOUTER,K),K=1,3)
 825        FORMAT (2(2X,F10.3,3X),2(2X,1PD10.3,3X),0PF11.4,
     .              3X,3(2X,1PD10.3,3X))
 827      CONTINUE
        END IF
C
C  back-substitution for all data arcs
        DO 850 I = 1 , NARCS
 830      READ   ( IFIN , 20 , END = 999 ) C133
          IF ( C133(42:47) .NE. 'UPDATE' ) GOTO 830
          call schrijf(iprint,ifout,c133)
 835      READ   ( IFIN , 20 , END = 999 ) C133
          IF ( C133(8:19) .NE. 'CROSS CORREL' ) THEN
            call schrijf(iprint,ifout,c133)
            GOTO 835
          END IF
 850    CONTINUE
      END IF
C-----------------------------------------------------------------------
C  station coordinates:
C  compute and print deltas during each outer iteration (except 1st)
      WRITE  ( IF6 , 860 )
 860  FORMAT ( /,'1CURRENT DELTAS AND STANDARD DEVIATIONS FOR ',
     .  'RECTANGULAR COORDINATES SOLUTIONS (CM):',
     .  //,'  station   iteration    delta(X)    delta(Y)    delta(Z)',
     .  '       sigma(X)    sigma(Y)    sigma(Z)',/,1X,100('-'),/)
C
      DO 900 IST = 1 ,NSTAT
        J = 1
        DO 870 K = 1 , 3
          DELTA(K) = ( XSTA(IST,J,K) - XAPSTA(IST,K) ) * 1.0D2
          SIGMA(K) = STDSTA(IST,J,K) * 1.0D2
 870    CONTINUE
        WRITE  ( IF6 , 880 ) ISTAT(IST) , J , ( DELTA(K) , K = 1 , 3 ) ,
     .    ( SIGMA(K) , K = 1 , 3  )
 880    FORMAT ( / , 1I7 , 1I10 , 2X , 3F12.2 , 3X , 3F12.2 )
C
        DO 895 J = 2 , IOUTER
          DO 885 K = 1 , 3
            DELTA(K) = ( XSTA(IST,J,K) - XSTA(IST,J-1,K) ) * 1.0D2
            SIGMA(K) = STDSTA(IST,J,K) * 1.0D2
 885      CONTINUE
          WRITE  ( IF6 , 890 ) J , ( DELTA(K) , K = 1 , 3 ) ,
     .      ( SIGMA(K) , K = 1 , 3 )
 890      FORMAT ( 7X , 1I10 , 2X , 3F12.2 , 3X , 3F12.2 )
 895    CONTINUE
 900  CONTINUE
C
C  earth rotation parameters:
C  compute and print deltas during each outer iteration (except 1st)
      IF ( NERP .NE. 0 ) THEN
        WRITE  ( IF6 , 910 )
 910    FORMAT ( /,'1CURRENT DELTAS AND STANDARD DEVIATIONS FOR ' ,
     .    'EARTH ROTATION PARAMETER SOLUTIONS ' ,
     .    '(MARCSEC, MARCSEC, MSEC):' ,
     .    // , '     begin       end         iteration   ' ,
     .    'delta(X)    delta(Y)  delta(A1-UT1)    ' ,
     .    'sigma(X)    sigma(Y)  sigma(A1-UT1)' ,
     .    / , 1X , 115('-') , / )
      END IF
C
      DO 950 IERP = 1 , NERP
        J = 1
        DO 920 K = 1 , 3
          DELTA(K) = ( XERP(IERP,J,K) - XAPERP(IERP,K) ) * 1.0D3
          SIGMA(K) = STDERP(IERP,J,K) * 1.0D3
 920    CONTINUE
        WRITE  ( IF6 , 930 ) T1ERP(IERP) , T2ERP(IERP) , J ,
     .    ( DELTA(K) , K = 1 , 3 ) , ( SIGMA(K) , K = 1 , 3 )
 930    FORMAT ( / , 2F12.3 , 1I10 , 2X , 3F12.4 , 3X , 3F12.4 )
C
        DO 945 J = 2 , IOUTER
          DO 935 K = 1 , 3
            DELTA(K) = ( XERP(IERP,J,K) - XERP(IERP,J-1,K) ) * 1.0D3
            SIGMA(K) = STDERP(IERP,J,K) * 1.0D3
 935      CONTINUE
          WRITE  ( IF6 , 940 ) J , ( DELTA(K) , K = 1 , 3 ) ,
     .      ( SIGMA(K) , K = 1 , 3 )
 940      FORMAT ( 24X , 1I10 , 2X , 3F12.4 , 3X , 3F12.4 )
 945    CONTINUE
 950  CONTINUE
      GOTO 9999
C
C998  WRITE  ( IF6 , * ) 'Unexpected end-of-file on SYSIN. Stop.'
C     GOTO 9999
 999  WRITE  ( IF6 , * ) 'Unexpected end-of-file on unit' ,
     .  IFIN , '. Stop.'
      GOTO 9999
C
  550 format (a)
 9999 END
C***********************************************************************
      SUBROUTINE SURGE ( IUNIT , C1 , C2 , * , * )
C
      CHARACTER * 1 CDUM
      CHARACTER * 4 C1,C2,CREAD
C
  10  READ   ( IUNIT , 20 , END = 40 ) CDUM , CREAD
  20  FORMAT ( 1A1 , 1A4 )
      IF ( CREAD .EQ. C2 ) RETURN 1
      IF ( CREAD .NE. C1 ) GOTO 10
      RETURN
C
  40  RETURN 2
C
      END
C***********************************************************************
      SUBROUTINE READI4 ( IF6 , C133 , I1 , I2 , IVALUE )
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
      CHARACTER * 133 C133
      CHARACTER * 300 FMAT
C
      DO 10 I = I1 , I2
        IF ( C133(I:I) .EQ. '*' ) THEN
          WRITE  ( IF6 , * )
          WRITE  ( IF6 , * )
     .      '* detected in field',I,' of the following line:'
          WRITE  ( IF6 , * ) C133
          WRITE  ( IF6 , * ) 'Dummy value 9999 assigned.'
          IVALUE = 9999
          RETURN
        END IF
  10  CONTINUE
C
      J1 = I1 - 1
      J2 = I2 - I1 + 1
      WRITE  ( FMAT , 20 ) J1 , J2
      READ   ( C133 , FMAT ) IVALUE
  20  FORMAT ( '(' , 1I3 , 'X , 1I' , 1I3 , ')' )
      RETURN
C
      END
C***********************************************************************
      SUBROUTINE READR8 ( IF6 , C133 , I1 , I2 , VALUE )
C
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )
      CHARACTER * 133 C133
      CHARACTER * 300 FMAT
C
      DO 10 I = I1 , I2
        IF ( C133(I:I) .EQ. '*' ) THEN
          WRITE ( IF6 , * )
          WRITE ( IF6 , * )
     .      '* detected in field',I,' of the following line:'
          WRITE ( IF6 , * ) C133
          WRITE ( IF6 , * ) 'Dummy value 9999.0 assigned.'
          VALUE = 9999.0D0
          RETURN
        END IF
  10  CONTINUE
C
      J1 = I1 - 1
      J2 = I2 - I1 + 1
      WRITE  ( FMAT , 20 ) J1 , J2
      READ   ( C133 , FMAT ) VALUE
  20  FORMAT ( '(' , 1I3 , 'X , 1F' , 1I3 , '.0 )' )
      RETURN
C
      END
C***********************************************************************
      subroutine schrijf(iprint,ifout,c133)
      character*133 c133
      if (iprint.ne.0) return
      l=lnblnk(c133)
      write (ifout,550) c133(:l)
  550 format (a)
      end
