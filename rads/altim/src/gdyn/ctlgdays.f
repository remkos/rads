**********************************************************************
*
*     This program can make an overview of:
*       - The total number of data passes per month per year,
*         taken at individual stations (you have to specify at
*         least the first date, first day of the year).
*       - The total number of data passes per week per time interval,
*         specifed by user, taken at individual stations. The weeks
*         run from Sunday, 0.0 hr, until the next Saturday, 24.0 hr.
*       - The total number of data passes per day per time interval,
*         specifed by user, taken at individual stations.
*
*     Input:
*     FT10  - The data catalogue (Kootwijk format) or measurements
*             in Geodyn Binary Format (GBF).
*     SYSIN - Input parameter
*                - 1=data catalogue en Kootwijk format
*                - 2=measurements in Geodyn Binary Format
*             Satellite id
*             Measurement type
*                - 20 = laser data
*                - 38 = DORIS data
*             Overview type
*                - 1 = monthly status
*                - 2 = weekly status
*                - 3 = daily status
*             Begin data in YYMMDD
*             End data in YYMMDD
*
*     (SYSIN data has to be on one record)
* 2  9105001  20  2 910101 921231
*     Output files:
*     none.
*
*     D.C. Kuijper, December 18, 1991.
*     R. Noomen, September 21, 1988.
*
**********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NST=50,NPS=80,NTOT=NST*NPS)
      DIMENSION NOPAS(NST,NPS),ISTAT(NST),CSTAT(NST)
      DIMENSION TIME1(NPS),TIME2(NPS)
      DATA ISTAT,NOPAS/NST*0,NTOT*0/
C
C  READ INPUT DATA
C
      READ (5,*) INPUT,IDSAT,IDTYP,ISELCT,IBEGDT,IENDDT
C
C  INITIAL SETTINGS
C
      NOYEAR=IBEGDT/10000
C
C  MONTHLY STATUS
C
      IF (ISELCT.EQ.1) THEN
        NOWEEK=12
        IF (IENDDT.NE.(IBEGDT+1130)) IENDDT=IBEGDT+1130
        DO 5 I=1,NOWEEK
          IYMD1=10000*NOYEAR+100*I+1
          IYMD2=10000*NOYEAR+100*I+101
          CALL MJDATE(2,MJD1,IYMD1,IFLUT1,IFLUT2,IFLUT3)
          CALL MJDATE(2,MJD2,IYMD2,IFLUT1,IFLUT2,IFLUT3)
          TIME1(I)=DBLE(MJD1)
          TIME2(I)=DBLE(MJD2)
  5     CONTINUE
        GOTO 25
      END IF
C
C  WEEKLY STATUS
C
      IF (ISELCT.EQ.2) THEN
        TWKREF=47240.0D0-7000.0D0
C  ALL WEEKS START ON SUNDAY 0.0 HR.
        TIMINT=7.0D0
        CALL MJDATE(2,MJDS,IBEGDT,IFLUT1,IFLUT2,IFLUT3)
        N=INT((DBLE(MJDS)-TWKREF)/TIMINT+1.0D-10)
        DMJDS=TWKREF+TIMINT*DBLE(N)-21.0D0
      END IF
C
C  DAYLY STATUS
C
      IF (ISELCT.EQ.3) THEN
        TIMINT=1.0D0
        IYMDS=IBEGDT-1
        CALL MJDATE(2,MJDSB,IYMDS,IFLUT1,IFLUT2,IFLUT3)
        CALL MJDATE(2,MJDSE,IENDDT,IFLUT1,IFLUT2,IFLUT3)
        IF ((MJDSE-MJDSB).GT.NPS) THEN
          MJDSE=MJDSB+NPS-1
          CALL MJDATE(1,MJDSE,IENDDT,IFLUT1,IFLUT2,IFLUT3)
        END IF
        DMJDS=1.D0*MJDSB
      END IF
C
      INDOLD=1
      INDNEW=1
      TIME1(INDOLD)=DMJDS
      TIME2(INDOLD)=DMJDS+TIMINT
C
  10  TIME1(INDNEW)=TIME1(INDOLD)+TIMINT
      TIME2(INDNEW)=TIME2(INDOLD)+TIMINT
      TIMEM=(TIME1(INDNEW)+TIME2(INDNEW))/2.0D0
      MJDM=TIMEM
      CALL MJDATE(1,MJDM,IYMDM,IYM,IMONTH,IFLUT3)
      IF (IYMDM.LT.IBEGDT) GOTO 10
      IF (IYMDM.GT.IENDDT) GOTO 20
C
      INDOLD=INDNEW
      INDNEW=INDNEW+1
      IF (INDNEW.GT.NPS) THEN
        WRITE(6,*)
     .  'THE NUMBER OF WEEKS OR DAYS EXCEEDS ',NPS,'. STOP.'
        STOP 99
      END IF
      GOTO 10
  20  NOWEEK=INDNEW-1
C
  25  CONTINUE
      IF (INPUT.EQ.1) THEN
        DO 30 I=1,NOWEEK
          TIME1(I)=TIME1(I)-33282.0D0
          TIME2(I)=TIME2(I)-33282.0D0
  30    CONTINUE
      END IF
C
      NOSTAT=0
      IF (INPUT.EQ.2)
     .  CALL READ2(IDTYP,IDSAT)
      END
C***********************************************************************
C
      SUBROUTINE READ2(IDTYP,IDSAT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NST=200,NPS=500)
      INTEGER*2 ITYPE
      DIMENSION IDATA(17),ITYPE(2)
      REAL*8 MJDBEG(NST,NPS),MJDEND(NST,NPS)
      DIMENSION ISTAT(NST),NOPAS(NST),NOOBS(NST,NPS)
      EQUIVALENCE (IDATA(1),ISAT),(IDATA(2),ITYPE(1)),
     .  (IDATA(3),ISTNEW),(IDATA(5),MJD),(IDATA(6),TFRAC)
      character*80 filenm
C
      DELTIM=60.0D0/1440.0D0
      NOSTAT=0
C
      call getarg(1,filenm)
      open (10,file=filenm,status='old',form='unformatted')
  10  READ (10,END=100) IDATA
      IF (ITYPE(1).NE.IDTYP) GOTO 10
      IF (IDSAT.NE.ISAT) GOTO 10
      DATNEW=MJD+TFRAC
C
  30  DO 40 I=1,NOSTAT
        IF (ISTNEW.EQ.ISTAT(I)) THEN
          N=NOPAS(I)
          IF (DABS(DATNEW-MJDBEG(I,N)).LT.DELTIM) THEN
            NOOBS(I,N)=NOOBS(I,N)+1
	    MJDEND(I,N)=DATNEW
            GOTO 10
          ENDIF
          N=N+1
      IF (NOSTAT.GT.NST) THEN
        WRITE (6,*) 'THE NUMBER OF PASSES EXCEEDS ',NPS,'. STOP.'
        STOP 99
      END IF
          NOOBS(I,N)=1
          MJDBEG(I,N)=DATNEW
          MJDEND(I,N)=DATNEW
          NOPAS(I)=N
          GOTO 10
        END IF
  40  CONTINUE
      NOSTAT=NOSTAT+1
      IF (NOSTAT.GT.NST) THEN
        WRITE (6,*) 'THE NUMBER OF STATIONS EXCEEDS ',NST,'. STOP.'
        STOP 99
      END IF
      I=NOSTAT
      ISTAT(I)=ISTNEW
      NOPAS(I)=1
      NOOBS(I,1)=1
      MJDBEG(I,1)=DATNEW
      MJDEND(I,1)=DATNEW
      GOTO 10
C
 100  CONTINUE
      DO 110 I=1,NOSTAT
         DO 110 N=1,NOPAS(I)
	    ISEC0=nint((mjdbeg(i,n)-46066)*86400)
	    ISEC1=nint((mjdend(i,n)-46066)*86400)
 110        WRITE (20,120) isec0,isec1,ISTAT(I),NOOBS(I,N)
 120  FORMAT (i10,i11,I5,I4)
*120  FORMAT (F11.5,F12.5,I5,I4)
      END
************************************************************************
C
      SUBROUTINE PRINST(INPUT,TIME1,TIME2,NOSTAT,ISTAT,NOWEEK,NOPAS,
     .  IDSAT,IDTYP,ISELCT,IBEGDT,IENDDT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NST=200,NPS=54)
      CHARACTER*2 C1,C2,CW
      CHARACTER*6 C6
      CHARACTER*12 FMT
      CHARACTER*130 C130,CHULP
      DIMENSION NOPAS(NST,NPS),ISTAT(NST),ISUM(NPS)
      DIMENSION TIME1(NPS),TIME2(NPS)
      DIMENSION IHUN(NPS),ITEN(NPS),IONE(NPS),C130(NST),IPR(NPS),
     .          C1(70),C2(0:100),C6(0:1000)
      DATA ISUM,NTOTOT/NPS*0,0/
C
      IF (ISELCT.EQ.1) THEN
        DO 10 I=0,999
  10      WRITE (C6(I),15) I
  15    FORMAT (I6)
      ELSE
        DO 20 I=0,99
  20      WRITE (C2(I),25) I
  25    FORMAT (I2)
      END IF
      DO 30 I=1,70
  30    C1(I)='--'
      C2(0)=' -'
      C2(100)='**'
      C6(0)='     -'
      C6(1000)='   ***'
C
      IF (ISELCT.EQ.3) THEN
        CALL MJDATE(4,IFLUT,IBEGDT,IYYB,IMMB,IDDB)
        CALL MJDATE(4,IFLUT,IENDDT,IYYE,IMME,IDDE)
        WRITE (6,1) IYYB,IMMB,IDDB
   1    FORMAT (/,' Begin Date: ',I2,'/',I2.2,'/',I2.2)
        WRITE (6,2) IYYE,IMME,IDDE
   2    FORMAT (' End Date  : ',I2,'/',I2.2,'/',I2.2)
      ELSE
        WRITE (6,3) IBEGDT/10000
   3    FORMAT (/,' Year      :     19',I2)
      END IF
      WRITE (6,4) IDSAT
   4  FORMAT (' Satellite :  ',I7)
      WRITE (6,5) IDTYP
   5  FORMAT (' Type      :  ',I7)
C
      IF (ISELCT.EQ.1) THEN
        WRITE (6,35)
  35    FORMAT(//,' Station      Jan   Feb   Mar   Apr   May   Jun',
     .  '   Jul   Aug   Sep   Oct   Nov   Dec',/,1X,130('-'),/)
        IND=12
        GOTO 55
      ELSE
        WRITE (6,40)
  40    FORMAT (//,' Data arc information:',//,
     .   '   Arc           Start         Stop',/,1X,37('-'),/)
        IND=53
      END IF
      DO 42 I=1,NOWEEK
        IF (INPUT.EQ.1) THEN
          MJD1=TIME1(I)+33282.0D0
          MJD2=TIME2(I)+33282.0D0
        ELSE
          MJD1=TIME1(I)
          MJD2=TIME2(I)
        END IF
        CALL MJDATE(1,MJD1,IYMD1,IFLUT1,IFLUT2,IFLUT3)
        CALL MJDATE(1,MJD2,IYMD2,IFLUT1,IFLUT2,IFLUT3)
        WRITE (6,41) I,IYMD1,IYMD2
  41    FORMAT (I5,10X,2(I6,'.00',5X))
  42  CONTINUE
C
      DO 50 I=1,NOWEEK
        ITEN(I)=DBLE(I)/10.0D0+1.0D-6
        IONE(I)=I-10*ITEN(I)
  50  CONTINUE
      WRITE (6,51) (ITEN(I),I=1,NOWEEK)
  51  FORMAT (//,' Station  ',99(I2))
      WRITE (6,52) (IONE(I),I=1,NOWEEK)
  52  FORMAT (10X,99(I2))
      WRITE (6,53) (C1(I),I=1,NOWEEK+8)
  53  FORMAT (1X,99(A2),/)
  55  CONTINUE
      DO 65 I=1,NOSTAT
        NTOTS=0
        DO 60 J=1,IND
          IPR(J)=NOPAS(I,J)
          IF (ISELCT.EQ.1.AND.NOPAS(I,J).GT.999) IPR(J)=1000
          IF (ISELCT.NE.1.AND.NOPAS(I,J).GT.99) IPR(J)=100
          NTOTS=NTOTS+NOPAS(I,J)
          NTOTOT=NTOTOT+NOPAS(I,J)
  60    CONTINUE
        IF (ISELCT.EQ.1) THEN
          WRITE (C130(I),61) (C6(IPR(J)),J=1,IND),NTOTS
  61      FORMAT (12(1A6),I8)
        ELSE
          WRITE(CW,'(I2)')NOWEEK
          FMT='('//CW//'(1A2),I7)'
          WRITE (C130(I),FMT) (C2(IPR(J)),J=1,NOWEEK),NTOTS
        END IF
        DO 64 J=1,IND
          ISUM(J)=ISUM(J)+NOPAS(I,J)
  64    CONTINUE
  65  CONTINUE
C
      DO 70 I=1,NOSTAT-1
      DO 70 J=I+1,NOSTAT
        IF (ISTAT(I).GT.ISTAT(J)) THEN
          IHULP=ISTAT(I)
          ISTAT(I)=ISTAT(J)
          ISTAT(J)=IHULP
          CHULP=C130(I)
          C130(I)=C130(J)
          C130(J)=CHULP
        END IF
  70  CONTINUE
      DO 73 I=1,NOSTAT
          WRITE (6,71) ISTAT(I),C130(I)
  71      FORMAT (I6,4X,1A112)
  73  CONTINUE
C
      IF (ISELCT.EQ.1) THEN
        WRITE (6,75) (ISUM(I),I=1,NOWEEK)
  75    FORMAT (/,'  Total : ',12(I6))
      ELSE
        DO 80 I=1,NOWEEK
          DHUN=DBLE(ISUM(I))/100.0D0+1.0D-6
          IHUN(I)=INT(DHUN)
          IHULP=ISUM(I)-100*IHUN(I)
          DTEN=DBLE(IHULP)/10.0D0+1.0D-6
          ITEN(I)=INT(DTEN)
          IONE(I)=IHULP-10*ITEN(I)
  80    CONTINUE
        WRITE (6,81) (IHUN(I),I=1,NOWEEK)
  81    FORMAT (/,'  Total : ',99(I2))
        WRITE (6,82) (ITEN(I),I=1,NOWEEK)
        WRITE (6,82) (IONE(I),I=1,NOWEEK)
  82    FORMAT (10X,99(I2))
      END IF
C
      WRITE (6,90) NTOTOT
  90  FORMAT (//,' Total number of passes: ',I8)
C
      RETURN
      END
