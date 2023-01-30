      program ons2gbf

* This is an improved version of ONS2GBF.
* Improved means: more accessible. If requires not to make any links
* to in and output files and thus is more user friendly.
*
* Syntax: ons2gbf inputfile np.outputfile [ ql.outputfile ]
*
* If "ql.outputfile" it will be removed upon exit.
*
* 29-Jul-1998 - Remko Scharroo - Adopted from ons2gbf by Ron Noomen
*  4-Aug-1998 - Splitted MAXPAS into MAXPNP and MAXPQL
*-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER     OFFSET
      INTEGER * 2 J5,LEAP
      REAL        STD,TROP,IONO,DISPL,CMAS
      REAL TARRAY(2), DUM
      CHARACTER LINE*69
C
      LOGICAL HEADER,NEWPAS
C
      PARAMETER ( MAXSAT =   12 )
      PARAMETER ( N30    =   30 )
      PARAMETER ( MAXSTA =   50 )
      PARAMETER ( MAXPNP =20000 )
      PARAMETER ( MAXPQL = 1000 )
      PARAMETER ( N1000  = 1000 )
C
      DIMENSION NMEAS(8,9999),NTOT(8)
      DIMENSION XUTC(N1000),XGPSCO(N1000),TUTC(N1000)
      DIMENSION ISATNP(MAXPNP),ISTNP(MAXPNP)
      DIMENSION STRTNP(MAXPNP),ENDNP(MAXPNP),NOBSNP(MAXPNP)
      DIMENSION ISATQL(MAXPQL),ISTQL(MAXPQL)
      DIMENSION STRTQL(MAXPQL),ENDQL(MAXPQL),NOBSQL(MAXPQL)
      DIMENSION JRLEAP(N30),LLEAP(N30)
      DIMENSION ISTWAV(N1000),T1WAV(N1000),T2WAV(N1000),IWAV(N1000)
      DIMENSION IDATA(17),J5(2)
      DIMENSION ISTAT2(MAXSTA,MAXSAT),IPSNP2(MAXSTA,MAXSAT)
      DIMENSION IPSQL2(MAXSTA,MAXSAT)
      DIMENSION NOBNP2(MAXSTA,MAXSAT),NOBQL2(MAXSTA,MAXSAT)
C
C ** The dimension of the next two array's depends on    **
C ** the number of satellites (see also DATA-statement)  **
      DIMENSION ISAT(MAXSAT),CMASS(MAXSAT)
C **                                                     **
      DIMENSION ITWO(31),DMJDY(99)
      DIMENSION IL(67)
C
      EQUIVALENCE (IDATA(1),J1),   (IDATA(2),J5(1)),   (IDATA(3),J9),
     .  (IDATA(4),J13),   (IDATA(5),J17),    (IDATA(6),TFRAC),
     .  (IDATA(8),RANGE), (IDATA(10),J37),   (IDATA(11),J41),
     .  (IDATA(12),STD),  (IDATA(13),TROP),  (IDATA(14),J53),
     .  (IDATA(15),IONO), (IDATA(16),DISPL), (IDATA(17),CMAS)
      EQUIVALENCE (I25,K48), (I32,IQ25), (I37,IQ30), (I41,IQ34)
C
      COMMON/MEAS/NTOT,NMEAS
C
      DATA IFIN,IFTIME,IFLEAP,IFSYS/5,11,12,13/
      DATA ICUR,VLIGHT/1,299792458.0D0/
C
C ** The number of satellites (NSAT) must be equal to  **
C ** the parameter MAXSAT from the parameter list          **
C
C  satellite ids:
C    7603901      LAGEOS-1
C    9207002      LAGEOS-2
C    9105001      ERS-1
C    9205201      TOPEX/Poseidon
C    7501001      Starlette
C    8606101      Ajisai
C    6503201      
C    9400301      Meteor-3
C    9502101      ERS-2
C    9306102      Stella
C    9305401      GPS-35
C    9401601      GPS-36
      DATA ISAT/
     .     7603901,9207002,9105001,9205201,7501001,8606101,6503201,
     .     9400301,9502101 , 9306102 , 9305401 , 9401601 /
      DATA CMASS/
     .     0.251D0,0.251D0,0.000D0,0.000D0,0.075D0,1.010D0,0.000D0,
     .     0.000D0,0.000D0 , 0.075D0 , 0.0D0 , 0.0D0 /
      DATA IF6,IQLOUT,INPOUT/6,20,21/

* Check arguments

      call getarg(1,line)
      if (line(1:2).eq.'-h') then
         write (0,600)
600      format (
     |'ons2gbf - Convert SLR data in ONS format to GBF'//
     |'syntax: ons2gbf [np.outputfile [ql.outputfile]] <inputfile'//
     |'where'/
     |'inputfile    : input file of SLR data in ONS format'/
     |'  (input is read from standard input)'/
     |'np.outputfile: output file of SLR normal points in GBF format'/
     |'ql.outputfile: output file of SLR quick-look data in GBF format'/
     |'  (output files are optional: not generated when left blank)'//
     |'Further print output comes to standard output')
	 goto 9999
      endif
      call getarg(1,line)
      if (line.eq.' ') then
         open (inpout,status='scratch',form='unformatted' )
      else
         open (inpout,file=line,form='unformatted' )
      endif
      call getarg(2,line)
      if (line.eq.' ') then
         open (iqlout,status='scratch',form='unformatted' )
      else
         open (iqlout,file=line,form='unformatted' )
      endif
C
      NSAT = MAXSAT
      CALL INITAR(MAXSTA,MAXSAT,NMEAS,NTOT,ISTAT2,IPSNP2,IPSQL2,
     .            NOBNP2,NOBQL2)
C
      CALL INITDT(IFTIME,N1000,NUTC,XUTC,XGPSCO,TUTC)
C     CALL INITLP(IFLEAP,N30,NLEAP,JRLEAP,LLEAP)
      CALL INITWL(IFSYS,N1000,NWAV,ISTWAV,T1WAV,T2WAV,IWAV)
      CALL INITYR(DMJDY)
      CALL INIT2(31,ITWO)
C
      NPASNP=0
      NPASQL=0
      NREC = 0
      J5(1)=20
      J37=0
      J41=0
      IONO=0.0E0
      DISPL=0.0E0
      HEADER = .FALSE.
      NEWPAS = .FALSE.
      N1 = 0
      N2 = 0
      N3 = 0
C
C * CHECKING LINE-TYPES *
C   Ni = 1 : Line is data-type identifier
C        2 : Line is a headerline
C        3 : Line is a dataline
C        4 : Line contains a wrong character
C        5 : Checksum of the line incorrect
C * i=3 for the current line
C * If HEADER=.FALSE. the next data-block will be ignored
C
C * Reading lines from data-file *
  10  READ(IFIN,15,END=300) LINE
  15  FORMAT(A69)
      NREC = NREC + 1
      IF ( LINE .EQ. ' ' ) GOTO 10
      N1=N2
      N2=N3
C
C * Detection Normal Point or Quicklook data *
      IF (LINE(1:5).EQ.'88888') THEN
         IFOUT=IQLOUT
         N3=1
         HEADER = .TRUE.
         GOTO 10
      END IF
      IF (LINE(1:5).EQ.'99999') THEN
         IFOUT=INPOUT
         N3=1
         HEADER = .TRUE.
         GOTO 10
      END IF
      IF ( .NOT. HEADER ) GOTO 10
C
C * Character check *
      DO 20 N=1,69
         IF(LINE(N:N).GE.'0'.AND.LINE(N:N).LE.'9') GOTO 20
         IF(LINE(N:N).EQ.' ') GOTO 20
         IF(LINE(N:N).EQ.'-') GOTO 20
         N3=4
  20  CONTINUE
      IF (N3.EQ.4) THEN
         IF (N2.EQ.1) THEN
            WRITE(6,'(1X,''Data-block will be ignored, there is a '',
     .         ''wrong character in the header-line:''/1X,A69)') LINE
            HEADER=.FALSE.
            GOTO 10
         END IF
         IF (N2.GE.2) THEN
            WRITE(6,'(1X,''Data-line will be ignored, there is a '',
     .        ''wrong character in the line:''/1X,A69)') LINE
            GOTO 10
         END IF
      END IF
C
C * Detection Headerline *
      READ (LINE(1:7),30) KS
  30  FORMAT(1I7)
      DO 40 N=1,NSAT
         IF (KS.EQ.ISAT(N)) GOTO 70
  40  CONTINUE
      GOTO 200
C
C **     HEADER-LINE     **
C
  70  N3=2
      HEADER=.TRUE.
      NEWPAS =.TRUE.
      DO 71 N=1,69
         IF(LINE(N:N).EQ.'0') LINE(N:N) = ' '
  71  CONTINUE
C
C * Verification checksum header-line *
      READ(LINE(1:54),80) (IL(N),N=1,52),ISUM
  80  FORMAT(BZ,52I1,I2)
      ITOT=0
      DO 90 N=1,52
         ITOT=ITOT + IL(N)
  90  CONTINUE
      ITOT=MOD(ITOT,100)
      IF(ITOT.NE.ISUM .AND. ISUM.NE.0) THEN
         DO 92 N=1,54
            IF(LINE(N:N).EQ.' ') LINE(N:N) = '0'
  92     CONTINUE
         WRITE(6,'(1X,''The checksum is not correct in the header'',
     .        ''-line:'',1I3/1X,A69/1X,''Process continues'',
     .        '' and checksum must be:'',1I3,/1X)') ISUM,LINE,ITOT
         DO 95 N=1,54
            IF(LINE(N:N).EQ.'0') LINE(N:N) = ' '
  95     CONTINUE
      END IF
C
C * Reading data from headerline *
C
      READ (LINE(1:51),100) K1,K8,K10,K13,K17,K19,
     .  K21,K25,K33,K39,K43,K44,K45,K46,K47,K48
 100  FORMAT (BZ,1I7,1I2,1I3,1I4,1I2,1I2,
     .  1I4,1I8,1I6,1I4,5I1,1I4)
C
      GOTO 10
C
C **     DATA-LINE     **
C
 200  N3=3
C
C * Checking if pass has a correct header-line *
      IF(N2.LE.1) THEN
         WRITE (6,'(1X,''Data-block ignored because of '',
     .         ''missing header-line. Nrec: '',1I6)') NREC
         HEADER=.FALSE.
         GOTO 10
      END IF
      IF(N2.GE.4.AND.N1.LE.1) THEN
         WRITE (6,'(1X,''Data-block ignored because of '',
     .         ''error in header-line. Nrec: '',1I6)') NREC
         HEADER=.FALSE.
         GOTO 10
      END IF
      IF (N1.NE.1.AND.N2.EQ.2) THEN
         WRITE (6,'(1X,''Data-block ignored because of '',
     .         ''missing data-type identifier. Nrec: '',1I6)') NREC
         HEADER=.FALSE.
         GOTO 10
      END IF
      IF ( HEADER ) GOTO 205
      GOTO 10
C
C * Verification checksum data-line *
 205  IF(IFOUT.EQ.IQLOUT) THEN
         READ(LINE(1:69),210) (IL(N),N=1,67),ISUM
 210     FORMAT(67I1,I2)
         ITOT=0
         DO 220 N=1,67
            ITOT=ITOT + IL(N)
 220     CONTINUE
         GOTO 250
      ELSE
         READ(LINE(1:63),230) (IL(N),N=1,52),ISUM
 230     FORMAT(52I1,I2)
         ITOT=0
         DO 240 N=1,52
            ITOT=ITOT + IL(N)
 240     CONTINUE
      END IF
 250  ITOT=MOD(ITOT,100)
      IF(ITOT.NE.ISUM .AND. ISUM.NE.0) THEN
C        WRITE(6,'(1X,''Data line will be ignored, the checksum is '',
C    .        ''not correct in the following line:'',1I3,
C    .        /1X,A69/1X,''Process continues'',
C    .        '' and checksum must be:'',1I3,/1X)') ISUM,LINE,ITOT
C       GOTO 10
         WRITE(6,'(1X,''Warning: the checksum in the following line '',
     .        ''is not correct:'',1I3,
     .        /1X,A69/1X,''Process continues'',
     .        '' and checksum must be:'',1I3,/1X)') ISUM,LINE,ITOT
      END IF
C
C * Reading data from data-line *
      IF (IFOUT.EQ.INPOUT) THEN
         READ (LINE(1:47),265) I1A,I1B,I13A,I13B,I25,I32,I37,I41
  265    FORMAT (4I6,1I7,1I5,1I4,1I3)
      ELSE
         READ (LINE(1:47),270) I1A,I1B,I13A,I13B,IQ25,IQ30,IQ34,IQ37
  270    FORMAT (4I6,1I5,1I4,1I3,1I8)
         IF ( IQ37 .EQ. 0 ) IQ37 = K25
      END IF
C
C * Checking time scale *
      IF ( K44.NE.7 .AND. K44.NE.4 .AND. K44.NE.3 ) THEN
        NMEAS(1,K13)=NMEAS(1,K13)+1
        NTOT(1)=NTOT(1)+1
        GOTO 10
      END IF
C
C * Checking year of century *
      IF (K8.LT.57) THEN
        NMEAS(2,K13)=NMEAS(2,K13)+1
        NTOT(2)=NTOT(2)+1
        GOTO 10
      END IF
C
C * Adding 1 to number of NP-passes or QL-passes *
      IF(NEWPAS) THEN
        IF (IFOUT.EQ.IQLOUT) THEN
          CALL PLUS1 ( 6 , NPASQL , MAXPQL , 'MAIN  ' ,'NPASQL' )
        END IF
        IF (IFOUT.EQ.INPOUT.AND.K43.NE.0) THEN
          CALL PLUS1 ( 6 , NPASNP , MAXPNP , 'MAIN  ' ,'NPASNP' )
        END IF
        IF (IFOUT.EQ.INPOUT.AND.K43.EQ.0) THEN
          CALL PLUS1 ( 6 , NPASQL , MAXPQL , 'MAIN  ' ,'NPASQL' )
        END IF
      END IF
C
C * Calculate time in MJD *
      DATE=DBLE(I1A)/86400.0D1+DBLE(I1B)/86400.0D7
C  store date of first observation in block
      IF(NEWPAS) DATE1=DATE+DMJDY(K8)+DBLE(K10)
C
C * Check if 24 h. is crossed *
      IF(ABS(DATE+DMJDY(K8)+DBLE(K10)-DATE1).GT.0.5D0) K10 = K10 + 1
      DATEXX = DATE + DMJDY(K8) + DBLE ( K10 )
C
C * convert from UTC(BIPM) to UTC(USNO) *
      IF (K44.EQ.7) THEN
        DELTAT=DIFFT(DATEXX,ICUR,K13,N1000,NUTC,XUTC,XGPSCO,TUTC)
        DATE=DATE-DELTAT/86400.0D6
      END IF
C
C * convert from UTC(GPS) to UTC(USNO) *
      IF (K44.EQ.4) THEN
C       DO 280 J=1,NLEAP
C       IF(K8.GE.JRLEAP(J)) THEN
C         LEAP=LLEAP(J)
C       END IF
C280    CONTINUE
        DELTAT = DIFFT(DATEXX,ICUR,K13,N1000,NUTC,XUTC,XGPSCO,TUTC)
        CO     = COGPS(DATEXX,ICUR,K13,N1000,NUTC,XUTC,XGPSCO,TUTC)
C       DATE=DATE-LEAP/86400.0D0+CO/86400.0D6-DELTAT/86400.0D6
        DATE=DATE+CO/86400.0D6-DELTAT/86400.0D6
      END IF
      KK44=3
C
      J1=K1
C
C ** The parameter J5(2) must be equal to 23 at this moment **
C **               J5(2)= 10*2 + K44 ( ==> KK44 )           **
C **     The factor 2 in '10*2' indicates the time of       **
C **     laser firing (see also comment at program start)   **
C **     In MERIT-II format it used to be '10*1' to         **
C **     indicate the satellite time.                       **
C **     At this stage the factor KK44 must be equal to 3   **
C **     to indicate the time has been transformed to       **
C **     UTC(USNO).                                         **
C
      J5(2)=10*2+KK44
      J9=K13
      IF (I32.NE.0) THEN
        ITEST=5
      ELSE
C  No meteo data available; hence no tropospheric correction can be
C  applied. Contrary to the description of the ONSITE data format,
C  the observation is already corrected for this delay.
        ITEST=0
      END IF
      IF (K1.EQ.9105001.OR.K1.EQ.9205201.OR.K1.EQ.9400301.OR.
     .    K1.EQ.9502101.OR.K1.EQ.9305401.OR.K1.EQ.9401601) THEN
        OFFSET=1
      ELSE
        OFFSET=0
      END IF
      J13=OFFSET*ITWO(28)+3*ITWO(25)+ITEST*ITWO(19)+1*ITWO(18)
      J17=INT(DATEXX)
      JJ17=INT(DATE+1.0D0)
      TFRAC=DATE+1.0D0-DBLE(JJ17)
      RANGE=(DBLE(I13A)*1.0D6+DBLE(I13B))
      IF (IFOUT.EQ.IQLOUT) RANGE=RANGE-DBLE(IQ37)
      RANGE=RANGE*0.5D-12*VLIGHT
      DO 290 I=1,NSAT
      IF (K1.EQ.ISAT(I)) THEN
        RANGE=RANGE+CMASS(I)
        IF (K1.EQ.9105001.OR.K1.EQ.9205201.OR.K1.EQ.9400301.OR.
     .      K1.EQ.9502101) THEN
          CMAS=0.0D0
        ELSE
          CMAS=CMASS(I)
        END IF
        GOTO 295
      END IF
 290  CONTINUE
      NMEAS(7,K13)=NMEAS(7,K13)+1
      NTOT(7)=NTOT(7)+1
 295  STD=DBLE(I25)*0.5D-12*VLIGHT
      K21HLP = 0
      CALL WAVE(K21HLP,K13,DATEXX,IWVCUR,ISTCUR,T1CUR,T2CUR,
     .  N1000,NWAV,ISTWAV,T1WAV,T2WAV,IWAV)
      IF (K21.NE.0) THEN
C The next IF-loop changes the units of the wavelength from 0.1 to 1.0 nm
C depending on the data value as follows:
C value 3000 to 9999: units 0.1 nm; value 0001 to 2999: units 1.0 nm. 
        IF (K21 .LE. 2999) THEN
          K21=K21*10
        END IF
        IF ((NINT(DBLE(K21)*1.0D-1).NE.K21HLP).AND.NEWPAS) THEN
          WRITE(6,296) K13,NINT(DBLE(K21)*1.0D-1),K21HLP,NREC
 296      FORMAT(' WARNING: Wavelength on data record for ',I4,
     |' <',I4,'> and system.data <',I4,'> mismatch (rec: ',I5,')')
        END IF
        K21HLP=NINT(DBLE(K21)*1.0D-1)
      END IF
      TROP=(DBLE(K21HLP)*1.0D-3+99.0D0)*1.0D-30
      IH=I41
      IT=(I37+5)/10
      IP=(I32+5)/10
      J53=IH*ITWO(24)+IT*ITWO(12)+IP
C
C  * Write data to files and fill data-array for making catalog *
      IF (IFOUT.EQ.INPOUT.AND.K43.EQ.0) THEN
         WRITE (IQLOUT) IDATA
         ISATQL(NPASQL)=IDATA(1)
         ISTQL(NPASQL)=IDATA(3)
         IF(NEWPAS) THEN
           STRTQL(NPASQL)=DATEXX
           ENDQL(NPASQL)=DATEXX
         END IF
         IF(DATEXX.LT.STRTQL(NPASQL)) THEN
            STRTQL(NPASQL)=DATEXX
         END IF
         IF(DATEXX.GT.ENDQL(NPASQL)) THEN
            ENDQL(NPASQL)=DATEXX
         END IF
         NOBSQL(NPASQL)=NOBSQL(NPASQL)+1
      ELSE
         WRITE (IFOUT) IDATA
         IF(IFOUT.EQ.INPOUT) THEN
           ISATNP(NPASNP)=IDATA(1)
           ISTNP(NPASNP)=IDATA(3)
           IF(NEWPAS) THEN
             STRTNP(NPASNP)=DATEXX
             ENDNP(NPASNP)=DATEXX
           END IF
           IF(DATEXX.LT.STRTNP(NPASNP)) THEN
              STRTNP(NPASNP)=DATEXX
           END IF
           IF(DATEXX.GT.ENDNP(NPASNP)) THEN
              ENDNP(NPASNP)=DATEXX
           END IF
           NOBSNP(NPASNP)=NOBSNP(NPASNP)+1
         ELSE
           ISATQL(NPASQL)=IDATA(1)
           ISTQL(NPASQL)=IDATA(3)
           IF(NEWPAS) THEN
             STRTQL(NPASQL)=DATEXX
             ENDQL(NPASQL)=DATEXX
           END IF
           IF(DATEXX.LT.STRTQL(NPASQL)) THEN
              STRTQL(NPASQL)=DATEXX
           END IF
           IF(DATEXX.GT.ENDQL(NPASQL)) THEN
              ENDQL(NPASQL)=DATEXX
           END IF
           NOBSQL(NPASQL)=NOBSQL(NPASQL)+1
         END IF
      END IF
C
      NEWPAS =.FALSE.
C
      GOTO 10
C
 300  WRITE(6,'(//,1X,''TOTAL NUMBER OF LINES READ :'',I7)') NREC
      ENDFILE(IQLOUT)
      ENDFILE(INPOUT)
      CALL CATALO(NSAT,ISAT,MAXPNP,MAXPQL,MAXSTA,MAXSAT,
     .  NPASNP,ISATNP,ISTNP,STRTNP,ENDNP,NOBSNP,
     .  NPASQL,ISATQL,ISTQL,STRTQL,ENDQL,NOBSQL,
     .  ISTAT2,IPSQL2,IPSNP2,NOBQL2,NOBNP2)
      CALL PRINT
      WRITE ( IF6 , * )
      WRITE ( IF6 , * ) 'Normal end of program ons2gbf.'
C  for use on convex only:
C     DUM = ETIME ( TARRAY )
C     WRITE ( IF6, 801 ) TARRAY(1)
C801  FORMAT ( ' Elapsed CPU time is:', F8.3, ' seconds' )
C
9999  END
C*********************************************************************
C
      DOUBLE PRECISION FUNCTION DIFFT(DATE,ICUR,K13,
     .  N1000 , NUTC , XUTC , XGPSCO , TUTC )
C
C  This function computes the time difference between the UTC(USNO) and
C  the UTC(BIH) scale.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XUTC(N1000),XGPSCO(N1000),TUTC(N1000)
      COMMON/MEAS/NTOT(8),NMEAS(8,9999)
C
      IF (DATE.GE.TUTC(ICUR).AND.DATE.LE.TUTC(ICUR+1)) GOTO 100
      ICUR=ICUR+1
      IF (ICUR.NE.NUTC) THEN
        IF (DATE.GE.TUTC(ICUR).AND.DATE.LE.TUTC(ICUR+1)) GOTO 100
      END IF
      IF (DATE.LT.TUTC(1)) THEN
        NMEAS(3,K13)=NMEAS(3,K13)+1
        NTOT(3)=NTOT(3)+1
        DIFFT=XUTC(1)
        ICUR=1
        RETURN
      END IF
      IF (DATE.GT.TUTC(NUTC)) THEN
        NMEAS(4,K13)=NMEAS(4,K13)+1
        NTOT(4)=NTOT(4)+1
        DIFFT=XUTC(NUTC)
        ICUR=NUTC-1
        RETURN
      END IF
C
      DO 10 I=1,NUTC-1
      IF (DATE.GE.TUTC(I).AND.DATE.LE.TUTC(I+1)) THEN
        ICUR=I
        GOTO 100
      END IF
  10  CONTINUE
C
 100  DIFFT=XUTC(ICUR)+(XUTC(ICUR+1)-XUTC(ICUR))*(DATE-TUTC(ICUR))/
     .  (TUTC(ICUR+1)-TUTC(ICUR))
      IF (ICUR+1.GT.NUTC) ICUR=ICUR-1
      RETURN
C
      END
C*********************************************************************
      DOUBLE PRECISION FUNCTION COGPS(DATE,ICUR,K13,
     .  N1000,NUTC,XUTC,XGPSCO,TUTC)
C
C  This function computes the correction factor Co concerning the time
C  difference UTC(BIH)-UTC(GPS). UTC(BIH)-UTC(GPS)= -(leap sec.) + Co.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XUTC(N1000),XGPSCO(N1000),TUTC(N1000)
      COMMON/MEAS/NTOT(8),NMEAS(8,9999)
C
      IF (DATE.GE.TUTC(ICUR).AND.DATE.LE.TUTC(ICUR+1)) GOTO 100
      ICUR=ICUR+1
      IF (ICUR.NE.NUTC) THEN
        IF (DATE.GE.TUTC(ICUR).AND.DATE.LE.TUTC(ICUR+1)) GOTO 100
      END IF
      IF (DATE.LT.TUTC(1)) THEN
        NMEAS(5,K13)=NMEAS(5,K13)+1
        NTOT(5)=NTOT(5)+1
        COGPS=0.0D0
        ICUR=1
        RETURN
      END IF
      IF (DATE.GT.TUTC(NUTC)) THEN
        NMEAS(4,K13)=NMEAS(4,K13)-1
        NTOT(4)=NTOT(4)-1
        NMEAS(6,K13)=NMEAS(6,K13)+1
        NTOT(6)=NTOT(6)+1
        COGPS=XGPSCO(NUTC)
        ICUR=NUTC-1
        RETURN
      END IF
C
      DO 10 I=1,NUTC-1
      IF (DATE.GE.TUTC(I).AND.DATE.LE.TUTC(I+1)) THEN
        ICUR=I
        GOTO 100
      END IF
  10  CONTINUE
C
 100  COGPS=XGPSCO(ICUR)+(XGPSCO(ICUR+1)-XGPSCO(ICUR))*(DATE-
     .  TUTC(ICUR))/(TUTC(ICUR+1)-TUTC(ICUR))
      IF (ICUR+1.GT.NUTC) ICUR=ICUR-1
      IF (ABS(COGPS).LT.1.0D-10) THEN
        NMEAS(5,K13)=NMEAS(5,K13)+1
        NTOT(5)=NTOT(5)+1
      END IF
      RETURN
C
      END
C*********************************************************************
C
      SUBROUTINE INIT2(N,ITWO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION ITWO(N)
C
      ITWO(1)=2
      DO 10 I=2,N
      ITWO(I)=2*ITWO(I-1)
  10  CONTINUE
      RETURN
C
      END
C*********************************************************************
C
      SUBROUTINE INITDT(IFILE,N1000,NUTC,XUTC,XGPSCO,TUTC)
C
C  This subroutine fills three array's with the data from a file (IFILE)
C  These array's will be used later to compute the time differences
C  between the various time scales.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION I(2),X(3)
      DIMENSION XUTC(N1000),XGPSCO(N1000),TUTC(N1000)
      character*80 sysdata
C
      call checkenv('ALTIM',sysdata,l)
      sysdata(l+1:)='/data/tables/times.data'
      open (ifile,file=sysdata,status='old')

      NUTC=0
      DO 20 J=1,3
      READ (IFILE,10,END=999)
  10  FORMAT (1A80)
  20  CONTINUE
C
  30  READ (IFILE,*,END=100) I,X
      IF (DABS(X(1)).LT.1.0D-10.AND.DABS(X(2)).LT.1.0D-10) GOTO 30
      CALL PLUS1 ( 6 , NUTC , N1000 , 'INITDT' ,'NUTC  ' )
      TUTC(NUTC)=DBLE(I(2))
      XUTC(NUTC)=X(1)
      XGPSCO(NUTC)=X(3)
      GOTO 30
 100  CONTINUE
      close (ifile)
      RETURN
C
 999  WRITE (6,*) 'INCORRECT NUMBER OF HEADER RECORDS IN FILE WITH ',
     .  'UTC(BIH)-UTC(USNO) DATA. STOP.'
      STOP 99
C
      END
C*********************************************************************
C
      SUBROUTINE INITLP(IFILE,N30,NLEAP,JRLEAP,LLEAP)
C
C  This subroutine fills two array's with data from a data-file (IFILE)
C  These array's contain the leap-seconds and the corresponding year
C  of thos leap-seconds. This leap-second is the difference between
C  UTC(BIH) and UTC(GPS).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER LINE*75
C
      DIMENSION JRLEAP(N30),LLEAP(N30)
C
      NLEAP=0
  10  READ (IFILE,20,END=999) LINE
  20  FORMAT (1A75)
      IF (LINE(1:5).NE.'A1UTC') THEN
        GOTO 10
      ELSE
  30    NLEAP=NLEAP+1
        IF (NLEAP.GT.N30) THEN
          WRITE (6,*) 'THE NUMBER OF A1-UTC TIME-DIFFERENCES EXCEEDS',
     .      N30,'. STOP.'
          STOP 99
        END IF
  40    READ (IFILE,50,END=999) LINE
  50    FORMAT (1A75)
        IF (LINE(1:10).EQ.'           ') GOTO 40
        IF (LINE(1:5).EQ.'FLUXS') GOTO 100
        READ (LINE(1:14),60,END=100) JRLEAP(NLEAP),LLEAP(NLEAP)
  60    FORMAT(I2,10X,I2)
        LLEAP(NLEAP)=LLEAP(NLEAP)-19
        GOTO 30
 100    NLEAP=NLEAP-1
        RETURN
      END IF
C
 999  WRITE (6,*) 'INSUFFICIENT DATA IN FILE WITH ',
     .  'A1(USNO)-UTC(USNO) DATA. STOP.'
      STOP 99
C
      END
C
C*********************************************************************
      SUBROUTINE INITWL(IFILE,N1000,NWAV,ISTWAV,T1WAV,T2WAV,IWAV)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION ISTWAV(N1000),T1WAV(N1000),T2WAV(N1000),IWAV(N1000)
      character*80 sysdata
C
      call checkenv('ALTIM',sysdata,l)
      sysdata(l+1:)='/data/tables/system.data'
      open (ifile,file=sysdata,status='old')
      REWIND IFILE
C
      I=1
  10  READ (IFILE,20,END=100) IHULP1,IHULP2,IYMD1,IYMD2,IWAV(I)
  20  FORMAT (1I4,4X,1I1,31X,2(2X,1I6),1X,1I4)
      ISTWAV(I) = 10000 * IHULP2 + IHULP1
      CALL MJDATE(2,MJD1,IYMD1,IFLUT1,IFLUT2,IFLUT3)
      CALL MJDATE(2,MJD2,IYMD2,IFLUT1,IFLUT2,IFLUT3)
      T1WAV(I)=DBLE(MJD1)
      T2WAV(I)=DBLE(MJD2+1)
      CALL PLUS1 ( 6 , I , N1000 , 'INITWL' ,'ISTWAV' )
      GOTO 10
C
 100  NWAV=I-1
      close (ifile)
      RETURN
C
      END
C*********************************************************************
C
      SUBROUTINE INITYR(DMJDY)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION DMJDY(99)
C
      DO 10 I=57,99
      CALL MJDATE(2,MJD,10000*I+101,IFLUT1,IFLUT2,IFLUT3)
      DMJDY(I)=DBLE(MJD)-1.0D0
  10  CONTINUE
      RETURN
C
      END
C*********************************************************************
C
      SUBROUTINE PRINT
C
C  This subroutine prints a overview of the various errors concerning
C  the translated data.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/MEAS/NTOT(8),NMEAS(8,9999)
C
      IF (NTOT(1).NE.0) THEN
        WRITE (6,100)
 100    FORMAT (///,' Overview of measurements deleted because of ',
     .    'epoch time scale undefined in GEODYN format:',//,
     .    ' station     measurements',/,1X,24('-'),/)
        DO 130 I=1,9999
        IF (NMEAS(1,I).NE.0) WRITE (6,120) I,NMEAS(1,I)
 120    FORMAT (1I6,5X,1I10)
 130    CONTINUE
      END IF
C
      IF (NTOT(2).NE.0) THEN
        WRITE (6,220)
 220    FORMAT (///,' Overview of measurements deleted because of ',
     .    'incorrect year of century:',//,
     .    ' station     measurements',/,1X,24('-'),/)
        DO 230 I=1,9999
        IF (NMEAS(2,I).NE.0) WRITE (6,120) I,NMEAS(2,I)
 230    CONTINUE
      END IF
C
      IF (NTOT(3).NE.0) THEN
        WRITE (6,320)
 320    FORMAT (///,' WARNING: Overview of measurements corrected',
     .    ' for the difference between',/,
     .    ' UTC(BIH) and UTC(USNO),',
     .    ' dated before the first value available from the',/,
     .    ' data file (this value was used):',//,
     .    ' station     measurements',/,1X,24('-'),/)
        DO 330 I=1,9999
        IF (NMEAS(3,I).NE.0) WRITE (6,120) I,NMEAS(3,I)
 330    CONTINUE
      END IF
C
      IF (NTOT(4).NE.0) THEN
        WRITE (6,420)
 420    FORMAT (///,' WARNING: Overview of measurements corrected',
     .    ' for the difference between',/,
     .    ' UTC(BIH) and UTC(USNO),',
     .    ' dated after the final value available from the',/,
     .    ' data file (this value was used):',//,
     .    ' station     measurements',/,1X,24('-'),/)
        DO 430 I=1,9999
        IF (NMEAS(4,I).NE.0) WRITE (6,120) I,NMEAS(4,I)
 430    CONTINUE
      END IF
C
      IF (NTOT(5).NE.0) THEN
        WRITE (6,520)
 520    FORMAT (///,' WARNING: Overview of measurements corrected',
     .    ' for the difference between GPS',/,
     .    ' and UTC(BIH), dated before the first non-zero',
     .    ' value available',/,
     .    ' from the data file (Co=0 was used):',//,
     .    ' station     measurements',/,1X,24('-'),/)
        DO 530 I=1,9999
        IF (NMEAS(5,I).NE.0) WRITE (6,120) I,NMEAS(5,I)
 530    CONTINUE
      END IF
C
      IF (NTOT(6).NE.0) THEN
        WRITE (6,620)
 620    FORMAT (///,' WARNING: Overview of measurements corrected',
     .    ' for the difference between',/,
     .    ' UTC(BIH) and UTC(USNO)',
     .    ' and for the difference between the GPS',
     .    ' and UTC(BIH),',/,
     .    ' dated after the final value available',
     .    ' from the data file (final values',/,
     .    ' were used for Co and UTC(USNO)-UTC(BIH).):',//,
     .    ' station     measurements',/,1X,24('-'),/)
        DO 630 I=1,9999
        IF (NMEAS(6,I).NE.0) WRITE (6,120) I,NMEAS(6,I)
 630    CONTINUE
      END IF
C
      IF (NTOT(7).NE.0) THEN
        WRITE (6,720)
 720    FORMAT (///,' WARNING: Overview of measurements not yet',
     .    ' corrected for center-of-mass offset',
     .    ' (unknown satellite id.)',/,
     .    ' (GEODYN can NOT apply this correction!):',//,
     .    ' station     measurements',/,1X,24('-'),/)
        DO 730 I=1,9999
        IF (NMEAS(7,I).NE.0) WRITE (6,120) I,NMEAS(7,I)
 730    CONTINUE
      END IF
C
      IF (NTOT(8).NE.0) THEN
        WRITE (6,820)
 820    FORMAT (///,' WARNING: No proper wavelength could be',
     .    ' found for the following measurements',/,
     .    ' (a value of 600 nm was encoded in the translated',/,
     .    ' measurements):',//,
     .    ' station     measurements',/,1X,24('-'),/)
        DO 830 I=1,9999
        IF (NMEAS(8,I).NE.0) WRITE (6,120) I,NMEAS(8,I)
 830    CONTINUE
      END IF
C
      RETURN
C
      END
C*********************************************************************
C
      SUBROUTINE WAVE(INEW,ISTNEW,TNEW,IOLD,ISTOLD,T1OLD,T2OLD,
     .  N1000,NWAV,ISTWAV,T1WAV,T2WAV,IWAV)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MEAS/NTOT(8),NMEAS(8,9999)
      DIMENSION ISTWAV(N1000),T1WAV(N1000),T2WAV(N1000),IWAV(N1000)
C
      IF (ISTNEW.EQ.ISTOLD.AND.TNEW.GE.T1OLD.AND.TNEW.LE.T2OLD) THEN
        INEW=IOLD
        RETURN
      END IF
C
      DO 10 I=1,NWAV
      IF (ISTWAV(I).NE.ISTNEW) GOTO 10
      IF (TNEW.LE.T1WAV(I).OR.TNEW.GE.T2WAV(I)) GOTO 10
      INEW=IWAV(I)
      IOLD=IWAV(I)
      ISTOLD=ISTWAV(I)
      T1OLD=T1WAV(I)
      T2OLD=T2WAV(I)
      RETURN
  10  CONTINUE
C
      WRITE ( 6 , * ) 'WAVE: no entry could be found for station ' ,
     .                ISTNEW , ' for epoch ' , TNEW
      WRITE ( 6 , * ) 'default value (600) for wavelength assumed'
C
      INEW=600
      IOLD=0
      ISTOLD=0
      T1OLD=0.0D0
      T2OLD=0.0D0
      NMEAS(8,ISTNEW)=NMEAS(8,ISTNEW)+1
      NTOT(8)=NTOT(8)+1
      RETURN
C
      END
C***********************************************************************
C
      SUBROUTINE CATALO(NSAT,ISAT,MAXPNP,MAXPQL,MAXSTA,MAXSAT,
     .  NPASNP,ISATNP,ISTNP,STRTNP,ENDNP,NOBSNP,
     .  NPASQL,ISATQL,ISTQL,STRTQL,ENDQL,NOBSQL,
     .  ISTAT2,IPSQL2,IPSNP2,NOBQL2,NOBNP2)
C
C  This subroutine prints a catalog of the translated data.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL STNEW
      DIMENSION ISATNP(MAXPNP),ISTNP(MAXPNP)
      DIMENSION STRTNP(MAXPNP),ENDNP(MAXPNP),NOBSNP(MAXPNP)
      DIMENSION ISATQL(MAXPQL),ISTQL(MAXPQL)
      DIMENSION STRTQL(MAXPQL),ENDQL(MAXPQL),NOBSQL(MAXPQL)
      DIMENSION ISTAT2(MAXSTA,MAXSAT),IPSNP2(MAXSTA,MAXSAT)
      DIMENSION IPSQL2(MAXSTA,MAXSAT)
      DIMENSION NOBNP2(MAXSTA,MAXSAT),NOBQL2(MAXSTA,MAXSAT)
      DIMENSION ISAT(MAXSAT)
C
      CALL SORTPS(MAXPQL,NPASQL,ISATQL,ISTQL,STRTQL,ENDQL,NOBSQL)
      CALL SORTPS(MAXPNP,NPASNP,ISATNP,ISTNP,STRTNP,ENDNP,NOBSNP)
C
*     WRITE(6,'(/,1X,''TOTAL NUMBER OF PASSES :'',I6)') (NPASQL+NPASNP)
*     WRITE(6,'(3X,''SAT'',5X,''STAT'',3X,'' START-TIME'',7X,
*    .       ''STOP-TIME'',3X,''NO.OF'',2X,''DATA'')')
*     WRITE(6,'(18X,2(''YYMMDD HHMMSS'',3X),''OBS'',2X,''TYPE'',
*    .       /,1X,58(''-''))')
*     I=1
*     J=1
*  5  IF(I.GT.NPASQL) THEN
*       IF(J.LE.NPASNP) THEN
*         CALL PRNTPS(ISATNP(J),ISTNP(J),STRTNP(J),ENDNP(J),
*    .          NOBSNP(J),'NP')
*         J=J+1
*       END IF
*     ELSE
*       IF(J.GT.NPASNP) THEN
*         CALL PRNTPS(ISATQL(I),ISTQL(I),STRTQL(I),ENDQL(I),
*    .          NOBSQL(I),'QL')
*         I=I+1
*       ELSE
*         IF(STRTQL(I).LT.STRTNP(J)) THEN
*           CALL PRNTPS(ISATQL(I),ISTQL(I),STRTQL(I),ENDQL(I),
*    .          NOBSQL(I),'QL')
*           I=I+1
*         ELSE
*           CALL PRNTPS(ISATNP(J),ISTNP(J),STRTNP(J),ENDNP(J),
*    .          NOBSNP(J),'NP')
*           J=J+1
*         END IF
*       END IF
*     END IF
*     IF(I.LT.(NPASQL+1).OR.J.LT.(NPASNP+1)) GOTO 5
*     WRITE(6,'(1X,58(''-''))')
C
C ** Block Quick-look data **
C
      WRITE(6,'(/,1X,''NUMBER OF QL PASSES :'',I5)') NPASQL
      WRITE(6,'(3X,''SAT'',5X,''STAT'',3X,'' START-TIME'',7X,
     .       ''STOP-TIME'',3X,''NO.OF'')')
      WRITE(6,'(18X,2(''YYMMDD HHMMSS'',3X),''OBS'',/,1X,52(''-''))')
      DO 60 I=1,NPASQL
        CALL PRNTPS(ISATQL(I),ISTQL(I),STRTQL(I),ENDQL(I),
     .          NOBSQL(I),'  ')
  60  CONTINUE
      WRITE(6,'(1X,52(''-''))')
C
C ** Block normal-points data **
      WRITE(6,'(/,1X,''NUMBER OF NP PASSES :'',I5)') NPASNP
      WRITE(6,'(3X,''SAT'',5X,''STAT'',3X,'' START-TIME'',7X,
     .       ''STOP-TIME'',3X,''NO.OF'')')
      WRITE(6,'(18X,2(''YYMMDD HHMMSS'',3X),''OBS'',/,1X,52(''-''))')
      DO 80 J=1,NPASNP
        CALL PRNTPS(ISATNP(J),ISTNP(J),STRTNP(J),ENDNP(J),
     .          NOBSNP(J),'  ')
  80  CONTINUE
      WRITE(6,'(1X,52(''-''))')
C
C *  PER STATION AND SATELLITE  *
      NSTAT = 0
      DO 130 J=1,NPASQL
        STNEW=.TRUE.
        DO 110 M=1,NSTAT
          DO 100 N=1,NSAT
            IF(ISTQL(J).EQ.ISTAT2(M,N)) THEN
              DO 90 K=1,NSAT
                IF(ISATQL(J).EQ.ISAT(K)) THEN
                  ISTAT2(M,K) = ISTQL(J)
                  IPSQL2(M,K) = IPSQL2(M,K) + 1
                  NOBQL2(M,K) = NOBQL2(M,K) + NOBSQL(J)
                  GOTO 130
                END IF
  90            CONTINUE
            END IF
 100        CONTINUE
 110      CONTINUE
        CALL PLUS1 ( 6 , NSTAT , MAXSTA , 'CATALO' ,'NSTAT' )
        DO 120 N=1,NSAT
          IF(ISATQL(J).EQ.ISAT(N)) THEN
            ISTAT2(NSTAT,N) = ISTQL(J)
            IPSQL2(NSTAT,N) = 1
            NOBQL2(NSTAT,N) = NOBSQL(J)
            GOTO 130
          END IF
 120      CONTINUE
 130    CONTINUE
C
      DO 155 I=1,NPASNP
        DO 150 M=1,NSTAT
          DO 140 N=1,NSAT
            IF(ISTNP(I).EQ.ISTAT2(M,N)) THEN
              DO 135 K=1,NSAT
                IF(ISATNP(I).EQ.ISAT(K)) THEN
                  ISTAT2(M,K) = ISTNP(I)
                  IPSNP2(M,K) = IPSNP2(M,K) + 1
                  NOBNP2(M,K) = NOBNP2(M,K) + NOBSNP(I)
                  GOTO 155
                END IF
 135            CONTINUE
            END IF
 140        CONTINUE
 150      CONTINUE
        CALL PLUS1 ( 6 , NSTAT , MAXSTA , 'CATALO' ,'NSTAT' )
        DO 154 N=1,NSAT
          IF(ISATNP(I).EQ.ISAT(N)) THEN
            ISTAT2(NSTAT,N) = ISTNP(I)
            IPSNP2(NSTAT,N) = 1
            NOBNP2(NSTAT,N) = NOBSNP(I)
            GOTO 155
          END IF
 154      CONTINUE
 155    CONTINUE
      CALL SORTST(MAXSTA,MAXSAT,NSTAT,NSAT,ISTAT2,IPSQL2,
     .            IPSNP2,NOBQL2,NOBNP2)
C
      WRITE(6,160)
 160  FORMAT(//,1X,'PASS AND OBSERVATION COUNTS FOR DATATYPE: 20')
      WRITE(6,170)
 170  FORMAT(/,1X,'COUNTS PER STATION-SATELLITE COMBINATION.')
      WRITE(6,180)
 180  FORMAT(1X,'STATION',2X,'SATELITTE',7X,'NP',8X,'NP',
     |8X,'QL',8X,'QL')
      WRITE(6,190)
 190  FORMAT(24X,2('PASSES',5X,'OBS.',5X)/1X,58('-'))
      DO 310 M=1,NSTAT
       DO 300 N=1,NSAT
        IF(ISTAT2(M,N).NE.0) THEN
            WRITE(6,220) ISTAT2(M,N),ISAT(N),IPSNP2(M,N),NOBNP2(M,N),
     .                   IPSQL2(M,N),NOBQL2(M,N)
 220        FORMAT(I7,I11,4I10)
        END IF
 300   CONTINUE
 310  CONTINUE
      WRITE(6,'(1X,58(''-''))')
      WRITE(6,320)
 320  FORMAT(/,1X,'COUNTS PER STATION.')
      WRITE(6,330)
 330  FORMAT(2X,'Station',2X,'Passes',2X,'Observations',/,1X,30('-'))
      DO 360 M=1,NSTAT
        IPAS=0
        IOBS=0
        DO 340 N=1,NSAT
          IF(ISTAT2(M,N).NE.0) THEN
            IST=ISTAT2(M,N)
          END IF
          IPAS=IPAS+IPSNP2(M,N)+IPSQL2(M,N)
          IOBS=IOBS+NOBNP2(M,N)+NOBQL2(M,N)
 340    CONTINUE
        WRITE(6,350) IST,IPAS,IOBS
 350    FORMAT(I7,I9,I15)
 360  CONTINUE
      WRITE(6,370)
 370  FORMAT(1X,30('-'),//,1X,'COUNTS PER SATELLITE')
      WRITE(6,380)
 380  FORMAT(2X,'Satelitte',2X,'Passes',2X,'Observations',/,1X,32('-'))
      DO 410 N=1,NSAT
        IPAS=0
        IOBS=0
        DO 390 M=1,NSTAT
          IPAS=IPAS+IPSNP2(M,N)+IPSQL2(M,N)
          IOBS=IOBS+NOBNP2(M,N)+NOBQL2(M,N)
 390    CONTINUE
        IF(IPAS.NE.0) THEN
          WRITE(6,400) ISAT(N),IPAS,IOBS
 400      FORMAT(I10,I8,I15)
        END IF
 410  CONTINUE
      WRITE(6,420)
 420  FORMAT(1X,32('-'))
      IPAS=0
      IOBS=0
      DO 440 M=1,NSTAT
        DO 430 N=1,NSAT
          IPAS=IPAS+IPSNP2(M,N)+IPSQL2(M,N)
          IOBS=IOBS+NOBNP2(M,N)+NOBQL2(M,N)
 430    CONTINUE
 440  CONTINUE
      WRITE(6,450) IPAS
 450  FORMAT(/,1X,'Total Counts:',1I9,'  Passes')
      WRITE(6,460) IOBS
 460  FORMAT(1I23,'  Observations',/)
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE PRNTPS(ISAT,IST,STRT,END,NOBS,TP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER * 2 TP
C
      TIMES=STRT
      CALL TYMD(TIMES,IYMDS,DHMSS)
      TIMEE=END
      CALL TYMD(TIMEE,IYMDE,DHMSE)
      WRITE(6,10) ISAT,IST,IYMDS,INT(DHMSS),
     .          IYMDE,INT(DHMSE),NOBS,TP
  10      FORMAT(1X,I7,3X,I4,3X,I6,1X,I6,' - ',I6,1X,I6,3X,I3,'  ',A2)
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SORTPS(MAXPAS,NPAS,ISAT,IST,STRT,END,NOBS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ISAT(MAXPAS),IST(MAXPAS)
      DIMENSION STRT(MAXPAS),END(MAXPAS),NOBS(MAXPAS)
C
      DO 100 I = 2 , NPAS
        IND = I
        IF ( STRT(I) .LT. STRT(I-1) ) THEN
          STRTI = STRT(I)
 10       IND = IND - 1
          IF ( IND .EQ. 1 ) GOTO 20
          IF ( STRTI .LT. STRT(IND-1) ) GOTO 10
 20       ISATI = ISAT(I)
          ISTI  = IST(I)
          ENDI  = END(I)
          NOBSI = NOBS(I)
          DO 30 J = I , IND + 1 , -1
            ISAT(J) = ISAT(J-1)
            IST(J)  = IST(J-1)
            END(J)  = END(J-1)
            STRT(J) = STRT(J-1)
            NOBS(J) = NOBS(J-1)
 30         CONTINUE
          ISAT(IND) = ISATI
          IST(IND)  = ISTI
          END(IND)  = ENDI
          STRT(IND) = STRTI
          NOBS(IND) = NOBSI
        ENDIF
 100  CONTINUE

      RETURN

      END
C***********************************************************************
C
      SUBROUTINE SORTST(MAXSTA,MAXSAT,NSTAT,NSAT,ISTAT,IPSQL,
     .                  IPSNP,NOBQL,NOBNP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ISTAT(MAXSTA,MAXSAT),IPSQL(MAXSTA,MAXSAT)
      DIMENSION IPSNP(MAXSTA,MAXSAT)
      DIMENSION NOBQL(MAXSTA,MAXSAT),NOBNP(MAXSTA,MAXSAT)
C
      DO 10 N = 1 , NSAT
        IF ( ISTAT(1,N) .NE. 0 ) THEN
          ISPREV = ISTAT(1,N)
          GOTO 11
        ENDIF
  10    CONTINUE
  11  CONTINUE
      DO 100 I = 2 , NSTAT
        DO 20 N = 1 , NSAT
          IF ( ISTAT(I,N) .NE. 0 ) THEN
            ISCUR = ISTAT(I,N)
            GOTO 21
          END IF
  20      CONTINUE
  21    CONTINUE
        IND = I
        IF ( ISCUR .LT. ISPREV ) THEN
          ISCNEW = ISPREV
  30      DO 35 N = 1 , NSAT
            IST            = ISTAT(IND,N)
            ISTAT(IND,N)   = ISTAT(IND-1,N)
            ISTAT(IND-1,N) = IST
            IQL            = IPSQL(IND,N)
            IPSQL(IND,N)   = IPSQL(IND-1,N)
            IPSQL(IND-1,N) = IQL
            INP            = IPSNP(IND,N)
            IPSNP(IND,N)   = IPSNP(IND-1,N)
            IPSNP(IND-1,N) = INP
            NQL            = NOBQL(IND,N)
            NOBQL(IND,N)   = NOBQL(IND-1,N)
            NOBQL(IND-1,N) = NQL
            NNP            = NOBNP(IND,N)
            NOBNP(IND,N)   = NOBNP(IND-1,N)
            NOBNP(IND-1,N) = NNP
  35        CONTINUE
          IND = IND - 1
          IF ( IND .EQ. 1 ) GOTO 50
          DO 40 N = 1 , NSAT
            IF ( ISTAT(IND-1,N) .NE. 0 ) THEN
              ISPREV = ISTAT(IND-1,N)
              GOTO 41
            END IF
  40        CONTINUE
  41      CONTINUE
          IF ( ISCUR .LT. ISPREV ) GOTO 30
  50      ISCUR = ISCNEW
        ENDIF
        ISPREV = ISCUR
 100    CONTINUE
C
      RETURN
      END
C---*-$--1----*----2----*----3----*----4----*----5----*----6----*----7-<--*----8

      SUBROUTINE INITAR(MAXSTA,MAXSAT,NMEAS,NTOT,ISTAT2,IPSNP2,IPSQL2,
     .                  NOBNP2,NOBQL2)

C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION NMEAS(8,9999),NTOT(8)
      DIMENSION ISTAT2(MAXSTA,MAXSAT),IPSNP2(MAXSTA,MAXSAT)
      DIMENSION IPSQL2(MAXSTA,MAXSAT)
      DIMENSION NOBNP2(MAXSTA,MAXSAT),NOBQL2(MAXSTA,MAXSAT)
C
      DO 20 I = 1 , 8
        NTOT(I) = 0
        DO 10 J = 1 , 9999
 10       NMEAS(I,J) = 0
 20     CONTINUE

C
      DO 40 I = 1 , MAXSAT
        DO 30 J = 1 , MAXSTA
          ISTAT2(J,I) = 0
          IPSNP2(J,I) = 0
          IPSQL2(J,I) = 0
          NOBNP2(J,I) = 0
          NOBQL2(J,I) = 0
 30       CONTINUE
 40     CONTINUE
C
      END
C*********************************************************************          
C                                                                               
      SUBROUTINE MJDATE(IND, MJD, IDATE, IYY, IMM, IDD)                         
C                                                                               
C        THIS SUBROUTINE DOES CONVERSIONS BETWEEN DATE AND MJD                  
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION MONTH(13)                                                       
      DATA MONTH / 0,31,59,90,120,151,181,212,243,273,304,334,365 /             
      DATA YEAR, MJCENT, IH / 365.25D0, 15019, 100 /                            
C                                                                               
      GOTO (30, 10, 20, 10), IND                                                
C        CONVERT IDATE TO IYY, IMM, IDD                                         
  10  IMM = IDATE/IH                                                            
      IDD = IDATE - IH*IMM                                                      
      IYY = IMM/IH                                                              
      IMM = IMM - IH *IYY                                                       
      IF (IND.EQ.4) RETURN                                                      
C                                                                               
C        CONVERT IYY, IMM, IDD TO MJD                                           
  20  MJD = MJCENT + IDD + MONTH(IMM) + (1461*IYY - 1)/4                        
      IF (IYY.EQ.4*(IYY/4) .AND. IMM.GT.2) MJD = MJD + 1                        
      RETURN                                                                    
C                                                                               
C        CONVERT MJD TO IYY, IMM, IDD                                           
  30  XMJD = MJD - MJCENT                                                       
      IYY = INT(XMJD/YEAR+1.0D-10)                                              
      ID = INT(XMJD -DBLE(IYY)*YEAR+1.D0+1.0D-10)                               
      IMM = ID/30 + 1                                                           
      IF (IYY.EQ.4*(IYY/4) ) IF (ID-60) 40, 50, 70                              
  40  IDD = ID - MONTH(IMM)                                                     
      IF (IDD.GT.0) GOTO 80                                                     
  50  IMM = IMM - 1                                                             
      GOTO 40                                                                   
  70  ID = ID - 1                                                               
      GOTO 40                                                                   
C        CONVERT IYY, IMM, IDD TO IDATE                                         
  80  IDATE = (IYY*IH + IMM)*IH + IDD                                           
C                                                                               
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
      SUBROUTINE TYMD(T,IYMD,DHMS)                                              
C                                                                               
C  subroutine TYMD converts T (in mjd.dd) to IYMD (in yymmdd)                   
C  and DHMS (in hhmmss.ss)                                                      
C                                                                               
      IMPLICIT REAL * 8 (A-H,O-Z)                                               
C                                                                               
      MJD=T                                                                     
      CALL MJDATE(1,MJD,IYMD,IFLUT1,IFLUT2,IFLUT3)                              
      TFRAC=T-DFLOAT(MJD)                                                       
      IH=TFRAC*24.0D0                                                           
      IM=TFRAC*1440.0D0-DFLOAT(IH)*60.0D0                                       
      DSEC=TFRAC*86400.0D0-DFLOAT(IM)*60.0D0-DFLOAT(IH)*3600.0D0                
      DHMS=DFLOAT(IH)*10000.0D0+DFLOAT(IM)*100.0D0+DSEC                         
C                                                                               
      RETURN                                                                    
      END                                                                       
C================================================================               
      SUBROUTINE PLUS1 ( IFILE , N1 , N2 , CSUB , CPARM )                       
C                                                                               
C----------------------------------------------------------------               
C  PLUS1 increases the value of N1 with 1                                       
C----------------------------------------------------------------               
C                                                                               
      IMPLICIT DOUBLE PRECISION ( A - H , O - Z )                               
      CHARACTER * 6 CSUB , CPARM                                                
C                                                                               
      N1 = N1 + 1                                                               
      IF ( N1 .GT. N2 ) THEN                                                    
        WRITE  ( IFILE , 10 ) CSUB , CPARM , N2                                 
  10    FORMAT ( 1X , 1A6 , ': the number of values for parameter ' ,           
     .           1A6 , ' exceeds' , 1I5 , '. stop.' )                           
        STOP 99                                                                 
      END IF                                                                    
      RETURN                                                                    
C                                                                               
      END                                                                       
