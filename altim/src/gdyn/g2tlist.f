      PROGRAM G2TLIST
************************************************************************
* Read Latitude, Longitude and Height from Geodyn-2 Trajectory file    *
* and sends the output either to standard output (formatted) or makes  *
* an ODR file.                                                         *
*                                                                      *
* Remko Scharroo, Delft, 3 July 1991.                                  *
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION    BUFFER(2048)
      CHARACTER*80 DECK(250),SPEC*4,SATEL*8,arg
      INTEGER*4    HGT,REQDAT
      EQUIVALENCE  (BUFFER(49),DECK)
*
* Read G2T header buffer
*
      REWIND (30)
      READ (30,ERR=1300) BUFFER
      IF (BUFFER(1).NE.-9D9) GOTO 1310
      NA=NINT(BUFFER(2))
      IF (BUFFER(7).NE.1) GOTO 1340
      NWDSAT=NINT(BUFFER(8))
      NTIMBF=NINT(BUFFER(10))
      DTIME=BUFFER(19)
      IOFF=0
      DO 10 I=202,207
   10    IF (BUFFER(I).GT.0) IOFF=IOFF+1
      IF (BUFFER(208).LE.0 .OR. BUFFER(209).LE.0 .OR.
     .    BUFFER(210).LE.0) GOTO 1320
*
* Read Alphanumeric header(s)
*
      READ (30,END=1330) BUFFER
      IF (BUFFER(1).NE.-8D9) GOTO 1310
      DO 30 I=2,NA
   30    READ (30,END=1330)
*
* Process Data buffers
*
      reqdat=0
      reqtim=0
      call getarg(1,arg)
      if (arg.ne.' ') read (arg,*) reqdat
      call getarg(2,arg)
      if (arg.ne.' ') read (arg,*) reqtim
      TIME=0
      KREC=0
  100 READ (30,END=1330) BUFFER
      IF (BUFFER(1).EQ.9D9) GOTO 9999
*
* Store date and start-time
*
      IDATE=INT(BUFFER(2)/1D6)
      TIME=BUFFER(2)-IDATE*1D6
      IHOUR=INT(TIME/1D4)
      TIME=TIME-IHOUR*1D4
      IMIN=INT(TIME/1D2)
      TIME=TIME-IMIN*1D2
      TIME=IHOUR*3600+IMIN*60+TIME
      NTB=NINT(BUFFER(5))
      IBUF=2*NTIMBF+5
      DO 120 IDAT=1,NTB
      IF (REQDAT.EQ.0)
     .      WRITE (6,630) IDATE,TIME,
     .      (BUFFER(I),I=IBUF+IOFF+1,IBUF+IOFF+3),
     .      (BUFFER(I),I=IBUF+IOFF-5,IBUF+IOFF  )
      IF (REQDAT.EQ.IDATE .AND. ABS(TIME-REQTIM).LT.1D0)
     .      WRITE (6,630) IDATE,TIME,
     .      (BUFFER(I),I=IBUF+IOFF+1,IBUF+IOFF+3),
     .      (BUFFER(I),I=IBUF+IOFF-5,IBUF+IOFF  )
         IBUF=IBUF+NWDSAT
  120    TIME=TIME+DTIME
      GOTO 100
*
* Formats
*
  550 FORMAT (A)
  630 FORMAT (I6,F10.3,2F11.6,F12.3/'ELEMS110   300  ',3F20.6/
     .'ELEMS2',10X,3F20.9)
*
* Errors
*
 1300 WRITE (6,*) 'ERROR'
      GOTO 9999
*
 1310 WRITE (6,550) 'G2TLIST: FILE IS NOT A G2T FILE'
      GOTO 9999
*
 1320 WRITE (6,550)
     .'G2TLIST: LATITUDE, LONGITUDE OR HEIGHT NOT AVAILABLE'
      GOTO 9999
*
 1330 WRITE (6,550) 'G2TLIST: PREMATURE END OF FILE'
      GOTO 9999
*
 1340 WRITE (6,550) 'G2TLIST: MULTIPLE SATELLITES NOT ALLOWED'
      GOTO 9999
*
 1350 WRITE (6,550) 'USAGE: G2TLIST <G2T-FILENAME>'//
     .' [ <ODR-FILENAME> <ARCNR> <REPEAT> <SATELLITE> ]'
 9999 END
