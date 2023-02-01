c irisub.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
C
C Includes subroutines IRISUB and IRIWEB to compute IRI parameters
C for specified location, date, time, and altitude range and subroutine
C IRIWEB to computes IRI parameters for specified location, date, time
C and variable range; variable can be altitude, latitude, longitude,
C year, month, day of month, day of year, or hour (UT or LT).
C IRIWEB requires IRISUB. Both subroutines require linking with the
c following library files IRIFUN.FOR, IRITEC.FOR,
c CIRA.FOR, IGRF.FOR
c
c !!!!! Subroutine INITIALIZE has to be called once before calling
c !!!!! IRISUB. This is already included in subroutine IRIWEB which
c !!!!! calls IRISUB.
C
c i/o units:  6   messages (during execution) to monitor
c
c-----------------------------------------------------------------------
c
C CHANGES FROM  IRIS11.FOR  TO   IRIS12.FOR:
C    - CIRA-1986 INSTEAD OF CIRA-1972 FOR NEUTRAL TEMPERATURE
C    - 10/30/91 VNER FOR NIGHTTIME LAY-VERSION:  ABS(..)
C    - 10/30/91 XNE(..) IN CASE OF LAY-VERSION
C    - 10/30/91 CHANGE SSIN=F/T TO IIQU=0,1,2
C    - 10/30/91 Te > Ti > Tn ENFORCED IN FINAL PROFILE
C    - 10/30/91 SUB ALL NAMES WITH 6 OR MORE CHARACTERS
C    - 10/31/91 CORRECTED HF1 IN HST SEARCH:  NE(HF1)>NME
C    - 11/14/91 C1=0 IF NO F1-REGION
C    - 11/14/91 CORRECTED HHMIN AND HZ FOR LIN. APP.
C    -  1/28/92 RZ12=0 included
C    -  1/29/92 NEQV instead of NE between URSIF2 and URSIFO
C    -  5/ 1/92 CCIR and URSI input as in IRID12
C    -  9/ 2/94 Decimal month (ZMONTH) for IONCOM
C    -  9/ 2/94 Replace B0POL with B0_TAB; better annually
C    -  1/ 4/95 DY for h>hmF2
C    -  2/ 2/95 IG for foF2, topside; RZ for hmF2, B0_TAB, foF1, NmD
C    -  2/ 2/95 winter no longer exclusive for F1 occurrrence
C    -  2/ 2/95 RZ and IG incl as DATA statement; smooth annual var.
C CHANGES FROM  IRIS12.FOR  TO   IRIS13.FOR:
C    - 10/26/95 incl year as input and corrected MODA; nrm for zmonth
C    - 10/26/95 use TCON and month-month interpolation in foF2, hmF2
C    - 10/26/95 TCON only if date changes
C    - 11/25/95 take out logicals TOPSI, BOTTO, and BELOWE
C    - 12/ 1/95 UT_LT for (date-)correct UT<->LT conversion
C    - 12/22/95 Change ZETA cov term to cov < 180; use cov inst covsat
C    -  2/23/96 take covmax(R<150) for topside; lyear,.. for lt
C    -  3/26/96 topside: 94.5/BETA inst 94.45/..; cov -> covsat(<=188)
C    -  5/01/96 No longer DY for h>hmF2 (because of discontinuity)
C    - 12/01/96 IRIV13: HOUR for IVAR=1 (height)
C    -  4/25/97 D-region: XKK le 10 with D1 calc accordingly.
C    -  1/12/97 DS model for lower ion compoistion DY model
C    -  5/19/98 seamon=zmonth if lati>0; zmonth= ...(1.0*iday)/..
C    -  5/19/98 DY ion composition model below 300 km now DS model
C    -  5/19/98 DS model includes N+, Cl down to 75 km HNIA changed
C    -  5/28/98 User input for Rz12, foF1/NmF1, hmF1, foE/NmE, hmE
C    -  9/ 2/98 1 instead of 0 in MODA after UT_LT call
C    -  4/30/99 constants moved from DATA statement into program
C    -  4/30/99 changed konsol-unit to 13 (12 is for IG_RZ).
C    -  5/29/99 the limit for IG comp. from Rz12-input is 174 not 274
C    - 11/08/99 jf(18)=t simple UT to LT conversion, otherwise UT_LT
C CHANGES FROM  IRIS13.FOR  TO   IRISUB.FOR:
c-----------------------------------------------------------------------
C-Version-MM/DD/YY-Description (person reporting correction)
C 2000.01 05/09/00 B0_98 replaces B0_TAB and B1: 1.9/day to 2.6/night
C 2000.02 06/11/00 including new F1 and indermediate region
C 2000.03 10/15/00 include Scherliess-Fejer drift model
C 2000.04 10/29/00 include special option for D region models
C 2000.05 12/07/00 change name IRIS13 to IRISUB
C 2000.06 12/14/00 jf(30),outf(20,100),oarr(50)
C 2000.07 03/17/01 include Truhlik-Triskova Te model and IGRF
C 2000.08 05/07/01 include Fuller-Rowell-Condrescu storm model
C 2000.09 07/09/01 LATI instead of LAT1 in F00 call (M. Torkar)
C 2000.10 07/09/01 sdte instead of dte in ELTEIK call (P. Wilkinson)
C 2000.11 09/18/01 correct computation of foF2 for Rz12 user input
C 2000.12 09/19/01 Call APF only if different date and time (P. Webb)
c 2000.13 10/28/02 replace TAB/6 blanks, enforce 72/line (D. Simpson)
C 2000.14 11/08/02 change unit for message file to 11 (13 is Kp)
C 2000.15 01/27/03 change F1_prob output; Te-IK for fix h and ELTE(h)
C 2000.16 02/04/03 along<0 -> along=along+360; F1 occ for hmf1&foF1
C 2000.17 02/05/03 zyear =12.97 (Dec 31); idayy=#days per year
C 2000.18 02/06/03 jf(27) for IG12 user input; all F1 prob in oar
C 2000.19 07/14/04 covsat<188 instead of covsat=<f(IG)<188
C 2000.19 02/09/05 declare INVDIP as real                F. Morgan
C 2000.20 11/09/05 replace B0B1 with BCOEF               T. Gulyaeva
C 2005.01 11/09/05 new topside ion composition; F107D from file
C 2005.02 11/14/05 jf(18) now IGRF/POGO; dip,mlat now IGRF10;
C 2005.03 11/15/05 sunrise/sunset/night for D,E,F1,F2; UT_LT removed
C 2005.04 05/06/06 FIRI D-region option not tied to peak
C 2005.04 05/06/06 Spread-F included, NeQuick included
C 2005.05 01/15/07 NeQuick uses CCIR-M3000F2 even if user-hmF2
C 2007.00 05/18/07 Release of IRI-2007
C
C*****************************************************************
C********* INTERNATIONAL REFERENCE IONOSPHERE (IRI). *************
C*****************************************************************
C**************** ALL-IN-ONE SUBROUTINE  *************************
C*****************************************************************
C
       SUBROUTINE IRISUB(JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR)
C-----------------------------------------------------------------
C
C INPUT:  JMAG          =0 geographic   = 1 geomagnetic coordinates
C         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C         IYYYY         Year as YYYY, e.g. 1985
C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL
C                          HOURS
C         HEIBEG,       HEIGHT RANGE IN KM; maximal 100 heights, i.e.
C          HEIEND,HEISTP        int((heiend-heibeg)/heistp)+1.le.100
C
C    JF switches to turn off/on (.true./.false.) several options
C
C    i       .true.                  .false.          standard version
C    -----------------------------------------------------------------
C    1    Ne computed            Ne not computed                     t
C    2    Te, Ti computed        Te, Ti not computed                 t
C    3    Ne & Ni computed       Ni not computed                     t
C    4    B0 - Table option      B0 - Gulyaeva (1987)                t
C    5    foF2 - CCIR            foF2 - URSI                     false
C    6    Ni - DS-78 & DY-85     Ni - DS-95 & TTS-03             false
C    7    Ne - Tops: f10.7<188   f10.7 unlimited                     t
C    8    foF2 from model        foF2 or NmF2 - user input           t
C    9    hmF2 from model        hmF2 or M3000F2 - user input        t
C   10    Te - Standard          Te - Using Te/Ne correlation        t
C   11    Ne - Standard Profile  Ne - Lay-function formalism         t
C   12    Messages to unit 6     no messages                         t
C   13    foF1 from model        foF1 or NmF1 - user input           t
C   14    hmF1 from model        hmF1 - user input (only Lay version)t
C   15    foE  from model        foE or NmE - user input             t
C   16    hmE  from model        hmE - user input                    t
C   17    Rz12 from file         Rz12 - user input                   t
C   18    IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
C   19    F1 probability model   critical solar zenith angle (old)   t
C   20    standard F1            standard F1 plus L condition        t
C   21    ion drift computed     ion drift not computed          false
C   22    ion densities in %     ion densities in m-3                t
C   23    Te_tops (Aeros,ISIS)   Te_topside (Intercosmos)        false
C   24    D-region: IRI-95       Special: 3 D-region models          t
C   25    F107D from AP.DAT      F107D user input (oarr(41))         t
C   26    foF2 storm model       no storm updating                   t
C   27    IG12 from file         IG12 - user input					 t
C   28    spread-F probability 	 not computed                    false
C   29    IRI01-topside          new options as def. by JF(30)   false
C   30    IRI01-topside corr.    NeQuick topside model   	     false
C     (29,30) = (t,t) IRIold, (f,t) IRIcor, (f,f) NeQuick, (t,f) TTS
C   ------------------------------------------------------------------
C
C  Depending on the jf() settings additional INPUT parameters may
c  be required:
C
C       Setting              INPUT parameter
C    -----------------------------------------------------------------
C    jf(8)  =.false.     OARR(1)=user input for foF2/MHz or NmF2/m-3
C    jf(9)  =.false.     OARR(2)=user input for hmF2/km or M(3000)F2
C    jf(10 )=.false.     OARR(15),OARR(16)=user input for Ne(300km),
C       Ne(400km)/m-3. Use OARR()=-1 if one of these values is not
C       available. If jf(23)=.false. then Ne(300km), Ne(550km)/m-3.
C    jf(13) =.false.     OARR(3)=user input for foF1/MHz or NmF1/m-3
C    jf(14) =.false.     OARR(4)=user input for hmF1/km
C    jf(15) =.false.     OARR(5)=user input for foE/MHz or NmE/m-3
C    jf(16) =.false.     OARR(6)=user input for hmE/km
C    jf(17) =.false.     OARR(33)=user input for Rz12
C    jf(21) =.true.      OARR(41)=user input for daily F10.7 index
C    jf(23) =.false.     OARR(41)=user input for daily F10.7 index
C    jf(24) =.false.     OARR(41)=user input for daily F10.7 index
C          optional for jf(21:24); default is F10.7D=COV
C    jf(25) =.false.     OARR(41)=user input for daily F10.7 index
C          if oarr(41).le.0 then 12-month running mean is
C          taken from internal file]
C    jf(27) =.false.     OARR(39)=user input for IG12
C    jf(21) =.true.      OARR(41)=user input for daily F10.7 index
c-----------------------------------------------------------------------
C*****************************************************************
C*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
C***     ELECTRON DENSITY         60/80 KM       1000 KM       ***
C***     TEMPERATURES              120 KM        2500/3000 KM  ***
C***     ION DENSITIES             100 KM        1000 KM       ***
C*****************************************************************
C*****************************************************************
C*********            INTERNALLY                    **************
C*********       ALL ANGLES ARE IN DEGREE           **************
C*********       ALL DENSITIES ARE IN M-3           **************
C*********       ALL ALTITUDES ARE IN KM            **************
C*********     ALL TEMPERATURES ARE IN KELVIN       **************
C*********     ALL TIMES ARE IN DECIMAL HOURS       **************
C*****************************************************************
C*****************************************************************
C*****************************************************************
      INTEGER	DAYNR
      REAL	LATI,LONGI,MODIP,NMF2,MAGBR,
     &		NMF1,NME,NMD,MLAT,MLONG
      CHARACTER	FILNAM*80,DIRNAM*80

      DIMENSION	ARIG(3),RZAR(3),F(3),E(4),XDELS(4),DNDS(4),
     &  FF0(988),XM0(441),F2(13,76,2),FM3(9,49,2),
     &  FF0N(988),XM0N(441),F2N(13,76,2),FM3N(9,49,2),INDAP(13)

      LOGICAL	EXT,SCHALT,sam_mon,sam_yea,sam_ut,
     &  sam_date,F1REG,sam_doy,dnight,enight,fnight

      COMMON	/CONST/UMR
      common	/ARGEXP/ARGMAX
      common	/BLOCK1/HMF2,NMF2,HMF1,F1REG
      common	/BLOCK2/B0,B1,C1
      common	/BLOCK3/HZ,T,HST
      common	/BLOCK4/HME,NME,HEF
      common	/BLOCK5/ENIGHT,E
      common	/BLOCK6/HMD,NMD,HDX
      common	/BLOCK7/D1,XKK,FP30,FP3U,FP1,FP2
      common	/BLO10/BETA,ETA,DELTA,ZETA
      common	/BLO11/B2TOP,TC3,hcor1
      common	/iounit/konsol,kerror

      INTEGER	montho/-9999/,iyearo/-9999/,idaynro/-9999/
      integer	unit,freeunit,konsol/-1/,kerror/0/
c set konsol=-1 if you do not want the konsol information

      EXTERNAL	XE1,XE2,XE3,XE4,XE5,XE6

      data	argmax/88.0/,umr/1.74532925199e-2/
      data	xdels/3*5.,10./,dnds/.016,.01,2*.016/

      save

C FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
C AGNR=OUTPUT (OUTPUT IS DISPLAYED OR STORED IN FILE OUTPUT.IRI)...

C CALCULATION OF DAY OF YEAR OR MONTH/DAY AND DECIMAL YEAR
c
c  leap year rule: years evenly divisible by 4 are leap years, except
c  years also evenly divisible by 100 are not leap years, except years
c  also evenly divisible by 400 are leap years. The year 2000 is a 100
c  and 400 year exception and therefore it is a normal leap year.
c  The next 100 year exception will be in the year 2100!

      iyear=iyyyy
      if(iyear.lt.100) iyear=iyear+1900
      idayy=365
      if(iyear/4*4.eq.iyear) idayy=366    ! leap year

      if(MMDD.lt.0) then
         DAYNR=-MMDD
         call MODA(1,iyear,MONTH,IDAY,DAYNR,nrdaym)
      else
         MONTH=MMDD/100
         IDAY=MMDD-MONTH*100
         call MODA(0,iyear,MONTH,IDAY,DAYNR,nrdaym)
      endif

      ryear = iyear + (daynr-1.0)/idayy

C CALCULATION OF GEODETIC/GEOMAGNETIC COORDINATES (LATI, LONGI AND
C MLAT, MLONG), MAGNETIC INCLINATION (DIP), DIP LATITUDE (MAGBR)
C AND MODIFIED DIP (MODIP), ALL IN DEGREES

      IF(JMAG.GT.0) THEN
         MLAT=ALATI
         MLONG=ALONG
      ELSE
         LATI=ALATI
         LONGI=ALONG
      ENDIF
      CALL GEODIP(IYEAR,LATI,LONGI,MLAT,MLONG,JMAG)
      if (LONGI.lt.0.) LONGI = LONGI + 360. ! -180/180 to 0-360

      call igrf_dip(lati,longi,ryear,300.0,dip,magbr,modip)

      ABSLAT=ABS(LATI)
      ABSMLT=ABS(MLAT)
      ABSMDP=ABS(MODIP)
      ABSMBR=ABS(MAGBR)

C CALCULATION OF UT/LT  ...............

      IF(DHOUR.le.24.0) then
          HOUR=DHOUR				! dhour =< 24 is LT
          ut=hour-longi/15.
          if(ut.lt.0) ut=ut+24.
      else
          UT=DHOUR-25.			 	! dhour>24 is UT+25
          hour=ut+longi/15.	 		! hour becomes LT
          if(hour.gt.24.) hour=hour-24.
      endif

c season assumes equal length seasons (92 days) with spring
c (season=1) starting at day-of-year=47; for lati<0 adjustment
c for southern hemisphere is made.

      SEASON=(DAYNR+45.0)/92.0
      ISEASON=INT(SEASON)
      IF (LATI.LT.0.0) ISEASON=ISEASON-2
      IF (ISEASON.LT.1) ISEASON=ISEASON+4

C 12-month running mean sunspot number (rssn) and Ionospheric Global
C index (gind)

      sam_mon=(month.eq.montho)
      sam_yea=(iyear.eq.iyearo)
      sam_doy=(daynr.eq.idaynro)
      sam_date=(sam_yea.and.sam_doy)
      sam_ut=(ut.eq.ut0)
      if(sam_date) goto 2910

      call tcon(iyear,month,iday,daynr,rzar,arig,ttt,nmonth)
      if(nmonth.lt.0) return		! jump to end of program

      rssn=rzar(3)
      gind=arig(3)
      COV=63.75+RSSN*(0.728+RSSN*0.00089)
      COVSAT=cov
      if(covsat.gt.188.) covsat=188

C CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG), SUN DECLINATION ANGLE
C (SUNDEC),SOLAR ZENITH ANGLE AT NOON (XHINON) AND TIME OF LOCAL
C SUNRISE/SUNSET (SAX, SUX; dec. hours) AT 70 KM (D-REGION), 110 KM
C (E-REGION), 200 KM (F1-REGION), AND 500 KM (TE, TI).

2910  continue
      CALL SOCO(daynr,HOUR,LATI,LONGI,70.,SUNDEC,XHI1,SAX70,SUX70)
      CALL SOCO(daynr,HOUR,LATI,LONGI,110.,SUD1,XHI2,SAX110,SUX110)
      CALL SOCO(daynr,HOUR,LATI,LONGI,200.,SUD1,XHI,SAX200,SUX200)
      CALL SOCO(daynr,HOUR,LATI,LONGI,500.,SUD1,XHI3,SAX500,SUX500)
      CALL SOCO(daynr,12.0,LATI,LONGI,110.,SUNDE1,XHINON,SAX1,SUX1)
      DNIGHT=.FALSE.
      if(abs(sax70).gt.25.0) then
         if(sax70.lt.0.0) DNIGHT=.TRUE.
         goto 9334
      endif
      if(SAX70.le.SUX70) goto 1386
      if((hour.gt.sux70).and.(hour.lt.sax70)) dnight=.true.
      goto 9334
1386  IF((HOUR.GT.SUX70).OR.(HOUR.LT.SAX70)) DNIGHT=.TRUE.

9334  ENIGHT=.FALSE.
      if(abs(sax110).gt.25.0) then
         if(sax110.lt.0.0) ENIGHT=.TRUE.
         goto 8334
      endif
      if(SAX110.le.SUX110) goto 9386
      if((hour.gt.sux110).and.(hour.lt.sax110)) enight=.true.
      goto 8334
9386  IF((HOUR.GT.SUX110).OR.(HOUR.LT.SAX110)) ENIGHT=.TRUE.

8334  FNIGHT=.FALSE.
      if(abs(sax200).gt.25.0) then
         if(sax200.lt.0.0) FNIGHT=.TRUE.
         goto 1334
      endif
      if(SAX200.le.SUX200) goto 7386
      if((hour.gt.sux200).and.(hour.lt.sax200)) fnight=.true.
      goto 1334
7386  IF((HOUR.GT.SUX200).OR.(HOUR.LT.SAX200)) FNIGHT=.TRUE.

C CALCULATION OF ELECTRON DENSITY PARAMETERS................
C lower height boundary (HNEA), upper boundary (HNEE)

1334  continue
      HNEA=65.
      IF(DNIGHT) HNEA=80.
      HNEE=2000.

      DELA=4.32
      IF(ABSMDP.GE.18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)
      DELL=1+EXP(-(ABSLAT-20.)/10.)

c E peak critical frequency (foF2), density (NmE), and height (hmE)

      FOE=FOEEDI(COV,XHI,XHINON,ABSLAT)
      NME=1.24E10*FOE*FOE

      HME=110.0

c F2 peak critical frequency foF2, density NmF2, and height hmF2
c
C READ CCIR AND URSI COEFFICIENT SET FOR CHOSEN MONTH

      IF(sam_mon.AND.(nmonth.EQ.nmono).and.sam_yea) GOTO 4292
      unit=freeunit()
      dirnam='/user/altim'
      call checkenv('ALTIM',dirnam,l)
      IF(sam_mon) GOTO 4293

c the program expects the coefficients files in ASCII format; if you
C want to use the binary version of the coefficients, please use the
C the statements that are commented-out below and comment-out the
C ASCII-related statements.

      WRITE(FILNAM,104) dirnam(:l),'ccir',MONTH+10
104   FORMAT(a,'/data/iri2007/',a,i2.2,'.asc')
      if (konsol.ge.0) write (konsol,*) "Reading: ",filnam
      OPEN(unit,FILE=FILNAM,STATUS='OLD',ERR=8448)
      READ(unit,4689) F2,FM3
4689  FORMAT(1X,4E15.8)
      CLOSE(unit)

C then URSI if chosen ....................................

      WRITE(FILNAM,104) dirnam(:l),'ursi',MONTH+10
      if (konsol.ge.0) write (konsol,*) "Reading: ",filnam
      OPEN(unit,FILE=FILNAM,STATUS='OLD',ERR=8448)
      READ(unit,4689) F2
      CLOSE(unit)

C READ CCIR AND URSI COEFFICIENT SET FOR NMONTH, i.e. previous
c month if day is less than 15 and following month otherwise

4293  continue

c first CCIR ..............................................

      WRITE(FILNAM,104) dirnam(:l),'ccir',NMONTH+10
      if (konsol.ge.0) write (konsol,*) "Reading: ",filnam
      OPEN(unit,FILE=FILNAM,STATUS='OLD',ERR=8448)
      READ(unit,4689) F2N,FM3N
      CLOSE(unit)

C then URSI if chosen .....................................

      WRITE(FILNAM,104) dirnam(:l),'ursi',NMONTH+10
      if (konsol.ge.0) write (konsol,*) "Reading: ",filnam
      OPEN(unit,FILE=FILNAM,STATUS='OLD',ERR=8448)
      READ(unit,4689) F2N
      CLOSE(unit)

C LINEAR INTERPOLATION IN SOLAR ACTIVITY. IG12 used for foF2

      RR2=ARIG(1)/100.
      RR2N=ARIG(2)/100.
      RR1=1.-RR2
      RR1N=1.-RR2N
      DO I=1,76
         DO J=1,13
            K=J+13*(I-1)
            FF0N(K)=F2N(J,I,1)*RR1N+F2N(J,I,2)*RR2N
            FF0(K)=F2(J,I,1)*RR1+F2(J,I,2)*RR2
	 enddo
      enddo

      RR2=RZAR(1)/100.
      RR2N=RZAR(2)/100.
      RR1=1.-RR2
      RR1N=1.-RR2N
      DO I=1,49
         DO J=1,9
            K=J+9*(I-1)
            XM0N(K)=FM3N(J,I,1)*RR1N+FM3N(J,I,2)*RR2N
            XM0(K)=FM3(J,I,1)*RR1+FM3(J,I,2)*RR2
         enddo
      enddo

4292  fof2   =  FOUT(MODIP,LATI,LONGI,UT,FF0)
      fof2n  =  FOUT(MODIP,LATI,LONGI,UT,FF0N)
      xm3000 = XMOUT(MODIP,LATI,LONGI,UT,XM0)
      xm300n = XMOUT(MODIP,LATI,LONGI,UT,XM0N)
      midm=15
      if(month.eq.2) midm=14
      if (iday.lt.midm) then
         fof2  = fof2n + ttt * (fof2-fof2n)
         xm3000= xm300n+ ttt * (xm3000-xm300n)
      else
         fof2  = fof2  + ttt * (fof2n-fof2)
         xm3000= xm3000+ ttt * (xm300n-xm3000)
      endif

      HMF2=HMF2ED(MAGBR,RSSN,FOF2/FOE,XM3000)

c stormtime updating

      if(.not.sam_date.or..not.sam_ut)
     &   call apf(iyear,month,iday,ut,indap)
      if(indap(1).gt.-2) then
         call STORM(indap,lati,longi,1,cglat,int(ut),daynr,stormcorr)
         fof2=fof2*stormcorr
      endif
      NMF2=1.24E10*FOF2*FOF2
      nmono=nmonth
      MONTHO=MONTH
      iyearo=iyear
      idaynro=daynr
      ut0=ut

c topside profile parameters .............................

      COS2=COS(MLAT*UMR)
      COS2=COS2*COS2
      FLU=(COVSAT-40.0)/30.0
      EX=EXP(-MLAT/15.)
      EX1=EX+1
      EPIN=4.*EX/(EX1*EX1)
      ETA1=-0.02*EPIN
      ETA = 0.058798 + ETA1 -
     &   FLU * (0.014065  - 0.0069724 * COS2) +
     &   FOF2* (0.0024287 + 0.0042810 * COS2  - 0.0001528 * FOF2)

      ZETA = 0.078922 - 0.0046702 * COS2 -
     &  FLU * (0.019132  - 0.0076545 * COS2) +
     &  FOF2* (0.0032513 + 0.0060290 * COS2  - 0.00020872 * FOF2)

      BETA=-128.03 + 20.253 * COS2 -
     &  FLU * (8.0755  + 0.65896 * COS2) +
     &  FOF2* (0.44041 + 0.71458 * COS2 - 0.042966 * FOF2)

      Z=EXP(94.5/BETA)
      Z1=Z+1
      Z2=Z/(BETA*Z1*Z1)
      DELTA=(ETA/Z1-ZETA/2.0)/(ETA*Z2+ZETA/400.0)

c NeQuick topside parameters (use CCIR-M3000F2 even if user-hmF2)

      dNdHmx=-3.467+1.714*log(foF2)+2.02*log(xM3000)
      dNdHmx=exp(dNdHmx)*0.01
      B2bot=0.04774*FOF2*FOF2/dNdHmx
      b2k = 3.22-0.0538*foF2-0.00664*hmF2+0.113*hmF2/B2bot
     &	+0.00257*rssn
      ee=exp(2.0*(b2k-1.0))
      b2k=(b2k*ee+1.0)/(ee+1.0)
      B2TOP=b2k*B2bot

c F layer - bottomside thickness parameter B0 and shape parameters B1

      B1=HPOL(HOUR,1.9,2.6,SAX200,SUX200,1.,1.)
      B0 = B0_98(HOUR,SAX200,SUX200,SEASON,RSSN,LONGI,MODIP)

c F1 layer height hmF1, critical frequency foF1, peak density NmF1

      FOF1=FOF1ED(ABSMBR,RSSN,XHI)
c F1 layer thickness parameter c1
      c1 = f1_c1(modip,hour,sux200,sax200)
c F1 occurrence probability with Scotto et al. 1997
      call f1_prob(xhi,mlat,rssn,f1pbw,f1pbl)
c F1 occurrence probability Ducharme et al.
      f1reg=(f1pbw.ge.0.5)
      NMF1=1.24E10*FOF1*FOF1

c E-valley: depth, width, height of deepest point (HDEEP),
c height of valley top (HEF)

      XDEL=XDELS(ISEASON)/DELA
      DNDHBR=DNDS(ISEASON)/DELA
      HDEEP=HPOL(HOUR,10.5/DELA,28.,SAX110,SUX110,1.,1.)
      WIDTH=HPOL(HOUR,17.8/DELA,45.+22./DELA,SAX110,SUX110,1.,1.)
      DEPTH=HPOL(HOUR,XDEL,81.,SAX110,SUX110,1.,1.)
      DLNDH=HPOL(HOUR,DNDHBR,.06,SAX110,SUX110,1.,1.)
      IF(DEPTH.GE.1.0) THEN
         IF(ENIGHT) DEPTH=-DEPTH
         CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
         IF(.NOT.EXT) GOTO 667
         if (konsol.ge.0) WRITE(KONSOL,650)
650      FORMAT(1X,'*NE* E-REGION VALLEY CAN NOT BE MODELLED')
      ENDIF
      WIDTH=.0
667   HEF=HME+WIDTH
      hefold=hef
      VNER = (1. - ABS(DEPTH) / 100.) * NME

c Parameters below E  .............................

      hmex=hme-9.
      NMD=XMDED(XHI,RSSN,4.0E8)
      HMD=HPOL(HOUR,81.0,88.0,SAX70,SUX70,1.,1.)
      F(1)=HPOL(HOUR,0.02+0.03/DELA,0.05,SAX70,SUX70,1.,1.)
      F(2)=HPOL(HOUR,4.6,4.5,SAX70,SUX70,1.,1.)
      F(3)=HPOL(HOUR,-11.5,-4.0,SAX70,SUX70,1.,1.)
      FP1=F(1)
      FP2=-FP1*FP1/2.0
      FP30=(-F(2)*FP2-FP1+1.0/F(2))/(F(2)*F(2))
      FP3U=(-F(3)*FP2-FP1-1.0/F(3))/(F(3)*F(3))
      HDX=HMD+F(2)

c indermediate region between D and E region; parameters xkk
c and d1 are found such that the function reaches hdx/xdx/dxdh

      X=HDX-HMD
      XDX=NMD*EXP(X*(FP1+X*(FP2+X*FP30)))
      DXDX=XDX*(FP1+X*(2.0*FP2+X*3.0*FP30))
      X=HME-HDX
      XKK=-DXDX*X/(XDX*ALOG(XDX/NME))

c if exponent xkk is larger than xkkmax, then xkk will be set to
c xkkmax and d1 will be determined such that the point hdx/xdx is
c reached; derivative is no longer continuous.

      xkkmax=5.
      if(xkk.gt.xkkmax) then
         xkk=xkkmax
         d1=-alog(xdx/nme)/(x**xkk)
      else
         D1=DXDX/(XDX*XKK*X**(XKK-1.0))
      endif

C SEARCH FOR HMF1 ..................................................

      hmf1=0
      IF(.not.F1REG) GOTO 380

c omit F1 feature if nmf1*0.9 is smaller than nme
      bnmf1=0.9*nmf1
      if(nme.ge.bnmf1) goto 9427

9245  XE2H=XE2(HEF)
      if(xe2h.gt.bnmf1) then
         hef=hef-1
         if(hef.le.hme) goto 9427
         goto 9245
      endif
      CALL REGFA1(HEF,HMF2,XE2H,NMF2,0.001,NMF1,XE2,SCHALT,HMF1)
      IF(.not.SCHALT) GOTO 3801

c omit F1 feature ....................................................

9427  if (konsol.ge.0) WRITE(KONSOL,11)
11    FORMAT(1X,'*NE* HMF1 IS NOT EVALUATED BY THE FUNCTION XE2'/
     &      1X,'CORR.: NO F1 REGION, B1=3, C1=0.0')
      HMF1=0.
      F1REG=.FALSE.

c Determine E-valley parameters if HEF was changed

3801  continue
      if(hef.ne.hefold) then
         width=hef-hme
         IF(ENIGHT) DEPTH=-DEPTH
         CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
         IF(.NOT.EXT) GOTO 380
         if (konsol.ge.0) WRITE(KONSOL,650)
         WIDTH=.0
         hef=hme
         hefold=hef
         goto 9245
      endif

C SEARCH FOR HST [NE3(HST)=NME] ......................................

380   continue

      IF(F1REG) then
          hf1=hmf1
          xf1=nmf1
      else
          hf1=(hmf2+hef)/2.
          xf1=xe2(hf1)
      endif

      hf2=hef
      xf2=xe3(hf2)
      if(xf2.gt.nme) goto 3885

      CALL REGFA1(hf1,HF2,XF1,XF2,0.001,NME,XE3,SCHALT,HST)
      if(schalt) goto 3885

      HZ=(HST+HF1)/2.0
      D=HZ-HST
      T=D*D/(HZ-HEF-D)
      GOTO 4933

c assume linear interpolation between HZ and HEF ..................

3885  if (konsol.ge.0) WRITE(KONSOL,100)
100   FORMAT(1X,'*NE* HST IS NOT EVALUATED BY THE FUNCTION XE3')
      HZ=(HEF+HF1)/2.
      xnehz=xe3(hz)
      if (konsol.ge.0) WRITE(KONSOL,901) HZ,HEF
901   FORMAT(6X,'CORR.: LIN. APP. BETWEEN HZ=',F5.1,' AND HEF=',F5.1)
      T=(XNEHZ-NME)/(HZ-HEF)
      HST=-333.

4933  CONTINUE

C CALCULATION FOR THE REQUIRED HEIGHT RANGE.......................
C In the absence of an F1 layer hmf1=hz since hmf1 is used in XE

      IF(hmf1.le.0.0) HMF1=HZ
      RETURN

* Error returns

8448  WRITE(kerror,8449) FILNAM
8449  FORMAT(' Could not read file ',a)
      return
      END
