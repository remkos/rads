      SUBROUTINE IRISUB(JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR)
*
* "Light" version of subroutine IRIS13 from the original IRI95
* code.
*
* 28-Jun-2002 - Created by Remko Scharroo from IRI95 code
*
C INPUT:  JMAG          =0 geographic   = 1 geomagnetic coordinates
C         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C         IYYYY         Year as YYYY, e.g. 1985
C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL 
C                          HOURS
C-------------------------------------------------------------------
C*****************************************************************
C*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
C***     ELECTRON DENSITY         60/80 KM       1000 KM       ***
C***     TEMPERATURES              120 KM        3000 KM       ***
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
      INTEGER   DAYNR,SEASON
      REAL      LATI,LONGI,MODIP,NMF2,MAGBR
      REAL      NMF1,NME,NMD,MLAT,MLONG
      CHARACTER FILNAM*80,dirnam*80
      real*4	xdels(4)/5.,5.,5.,10./
      real*4	F(3),E(4),DNDS(4)/.016,.01,.016,.016/
      real*4	FF0(988),XM0(441),F2(13,76,2),FM3(9,49,2)
      real*4	FF0N(988),XM0N(441),F2N(13,76,2),FM3N(9,49,2)
      real*4	ARIG(3),RZAR(3)
      LOGICAL	EXT,SCHALT,NIGHT,sam_mon,sam_yea
      LOGICAL	F1REG
      LOGICAL	sam_doy
      COMMON    /BLOCK1/HMF2,NMF2,HMF1
      common	/CONST/UMR
      common	/const1/humr
      common    /BLOCK2/B0,B1,C1
      common	/BLOCK3/HZ,T,HST,STR
      common    /BLOCK4/HME,NME,HEF
      common	/BLOCK5/NIGHT,E
      common    /BLOCK6/HMD,NMD,HDX
      common	/BLOCK7/D1,XKK,FP30,FP3U,FP1,FP2
      common    /BLO10/BETA,ETA,DELTA,ZETA,EPTR1X0,EPTR2X0
      common    /ARGEXP/ARGMAX
      EXTERNAL  XE1,XE2,XE3,XE4,XE5,XE6
      real	montho/-9999/,iyearo/-9999/,idayno/-9999/
      integer	unit,freeunit,konsol/-1/
      save
C
C PROGRAM CONSTANTS
C
      argmax=88.0
      umr=atan(1.)/45.
      humr=umr*15.
C
C CALCULATION OF GEOG. OR GEOM. COORDINATES IN DEG....................
C CALCULATION OF MAGNETIC INCLINATION (DIP), DECLINATION (DEC)........
C   DIP LATITUDE (MAGBR) AND MODIFIED DIP (MODIP). ALL IN DEGREE......
C
      IF(JMAG.GT.0) THEN
         MLAT=ALATI
         MLONG=ALONG
      ELSE
         LATI=ALATI
         LONGI=ALONG
      ENDIF
      CALL GGM(JMAG,LONGI,LATI,MLONG,MLAT)
      ABSLAT=ABS(LATI)
      CALL FIELDG(LATI,LONGI,300.0,XMA,YMA,ZMA,BET,DIP,DEC,MODIP)
      MAGBR=ATAN(0.5*TAN(DIP*UMR))/UMR
      ABSMLT=ABS(MLAT)
      ABSMDP=ABS(MODIP)
      ABSMBR=ABS(MAGBR)
C
C CALCULATION OF DAY OF YEAR AND SUN DECLINATION......................
C CALCULATION OF UT/LT AND RELATED YEAR, MONTH, DAYNRs ...............
C CALCULATION OF (UT-)SEASON (SUMMER=2, WINTER=4).....................
C
      iyear=iyyyy
      if(iyear.lt.100) iyear=iyear+1900
      if(MMDD.lt.0) then
          DAYNR=-MMDD
          call MODA(1,iyear,MONTH,IDAY,DAYNR,nrdaym)
      else
          MONTH=MMDD/100
          IDAY=MMDD-MONTH*100
          call MODA(0,iyear,MONTH,IDAY,DAYNR,nrdaym)
      endif
c
c lyear,lmonth,lday,ldaynre,lnrday related to LT
c
      lyear=iyear
      lmonth=month
      lday=iday
      ldaynr=daynr
      lnrday=nrdaym

      if (dhour.gt.24.) then
         UT=DHOUR-25.
         iytmp=iyear
         idtmp=daynr
         hour=ut+longi/15.
         if(hour.gt.24.) hour=hour-24.
      else
         HOUR=DHOUR
         iytmp=lyear
         idtmp=ldaynr
         ut=hour-longi/15.
         if(ut.lt.0) ut=ut+24.
      endif

      SEASON=INT((DAYNR+45.0)/92.0)
      IF(SEASON.LT.1) SEASON=4
      NSEASN=SEASON
      IF(LATI.le.0.0) then
         SEASON=SEASON-2
         IF(SEASON.LT.1) SEASON=SEASON+4
      endif
C
C CALCULATION OF MEAN F10.7CM SOLAR RADIO FLUX (COV)................
C CALCULATION OF RESTRICTED SOLAR ACTIVITIES (RSAT,COVSAT)..............
C
      sam_mon=(month.eq.montho)
      sam_yea=(iyear.eq.iyearo)
      sam_doy=(daynr.eq.idayno)
      if (.not.sam_yea.or..not.sam_doy) then
         call tcon(iyear,month,iday,daynr,rzar,arig,ttt,nmonth)
         if(nmonth.lt.0) return
         rssn=rzar(3)
         rsat=arig(3)
         COV=63.75+RSSN*(0.728+RSSN*0.00089)
         rlimit=rsat
         COVSAT=63.75+rlimit*(0.728+rlimit*0.00089)
         if(covsat.gt.188.) covsat=188
      endif
C
C CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG).........................
C NOON VALUE (XHINON).................................................
C
      CALL SOCO(ldaynr,HOUR,LATI,LONGI,SUNDEC,XHI,SAX,SUX)
      CALL SOCO(ldaynr,12.0,LATI,LONGI,SUNDE1,XHINON,SAXNON,SUXNON)

      night=.false.
      if(abs(sax).gt.25.0) then
         if(sax.lt.0.0) night=.true.
      else if(sax.le.sux) then
         if(hour.gt.sux.or.hour.lt.sax) night=.true.
      else
         if(hour.gt.sux.and.hour.lt.sax) night=.true.
      endif
C
C CALCULATION OF ELECTRON DENSITY PARAMETERS................
C
      HNEA=65.
      IF(NIGHT) HNEA=80.
      HNEE=2000.
      DELA=4.32
      IF(ABSMDP.GE.18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)
      DELL=1+EXP(-(ABSLAT-20.)/10.)
C!!!!!!! F-REGION PARAMETERS AND E-PEAK !!!!!!!!!!!!!!!!!!!!!!!!!!
      FOE=FOEEDI(COV,XHI,XHINON,ABSLAT)
      NME=1.24E10*FOE*FOE
      HME=110.0
C
C READ CCIR AND URSI COEFFICIENT SET FOR CHOSEN MONTH ............
C
      IF(sam_mon.AND.(nmonth.EQ.nmono).and.sam_yea) GOTO 4292
      unit=freeunit()
      dirnam='/user/altim'
      call checkenv('ALTIM',dirnam,l)
      IF(sam_mon) GOTO 4293
c
c the program expects the coefficients files in ASCII format; if you
C want to use the binary version of the coefficients, please use the
C the statements that are commented-out below and comment-out the
C ASCII-related statements.
c
      if(konsol.ge.0)write(konsol,*)iyear,month,iday,filnam
      WRITE(FILNAM,104) dirnam(:l),'ccir',MONTH+10
      if(konsol.ge.0)write(konsol,*)filnam
104   FORMAT(a,'/data/iri95/',a4,i2,'.asc')
      OPEN(unit,FILE=FILNAM,STATUS='OLD',ERR=8448)
      READ(unit,4689) F2,FM3
4689  FORMAT(1X,4E15.8)
      CLOSE(unit)
C
C then URSI if chosen ....................................
C
      WRITE(FILNAM,104) dirnam(:l),'ursi',MONTH+10
      if(konsol.ge.0)write(konsol,*)filnam
      OPEN(unit,FILE=FILNAM,STATUS='OLD',ERR=8448)
      READ(unit,4689) F2
      CLOSE(unit)
C
C READ CCIR AND URSI COEFFICIENT SET FOR NMONTH, i.e. previous 
c month if day is less than 15 and following month otherwise 
C
4293  continue
c
c first CCIR ..............................................
c
      if(konsol.ge.0)write(konsol,*)iyear,month,iday
      WRITE(FILNAM,104) dirnam(:l),'ccir',NMONTH+10
      if(konsol.ge.0)write(konsol,*)filnam
      OPEN(unit,FILE=FILNAM,STATUS='OLD',ERR=8448)
      READ(unit,4689) F2N,FM3N
      CLOSE(unit)
C
C then URSI if chosen .....................................
C
      WRITE(FILNAM,104) dirnam(:l),'ursi',NMONTH+10
      if(konsol.ge.0)write(konsol,*)filnam
      OPEN(unit,FILE=FILNAM,STATUS='OLD',ERR=8448)
      READ(unit,4689) F2N
      CLOSE(unit)

      nmono=nmonth
      MONTHO=MONTH
      iyearo=iyear
      idayno=daynr
C
C LINEAR INTERPOLATION IN SOLAR ACTIVITY. RSAT used for foF2
C
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

4292  zfof2  = FOUT(MODIP,LATI,LONGI,UT,FF0)
      fof2n  = FOUT(MODIP,LATI,LONGI,UT,FF0N)
      zm3000 = XMOUT(MODIP,LATI,LONGI,UT,XM0)
      xm300n = XMOUT(MODIP,LATI,LONGI,UT,XM0N)

      midm=15
      if(month.eq.2) midm=14
      if (iday.lt.midm) then
         yfof2 = fof2n + ttt * (zfof2-fof2n)
         xm3000= xm300n+ ttt * (zm3000-xm300n)
      else
         yfof2 = zfof2 + ttt * (fof2n-zfof2)
         xm3000= zm3000+ ttt * (xm300n-zm3000)
      endif
      FOF2=YFOF2
      NMF2=1.24E10*FOF2*FOF2
      HMF2=HMF2ED(MAGBR,RSSN,FOF2/FOE,XM3000)
c
c topside profile parameters .............................
c
      COS2=COS(MLAT*UMR)
      COS2=COS2*COS2
      FLU=(COVSAT-40.0)/30.0
      EX=EXP(-MLAT/15.)
      EX1=EX+1
      EPIN=4.*EX/(EX1*EX1)
      ETA1=-0.02*EPIN
      ETA=0.058798+ETA1+FLU*(-0.014065+0.0069724*COS2)+
     &(0.0024287+0.0042810*COS2-0.00015280*FOF2)*FOF2
      ZETA=0.078922-0.0046702*COS2-0.019132*FLU+0.0076545*FLU*COS2+
     &(0.0032513+0.0060290*COS2-0.00020872*FOF2)*FOF2
      BETA=-128.03+20.253*COS2+FLU*(-8.0755-0.65896*COS2)+(0.44041
     &+0.71458*COS2-0.042966*FOF2)*FOF2
      Z=EXP(94.5/BETA)
      Z1=Z+1
      Z2=Z/(BETA*Z1*Z1)
      DELTA=(ETA/Z1-ZETA/2.0)/(ETA*Z2+ZETA/400.0)
      EPTR1X0=EPTR(300.-DELTA,BETA,394.5)
      EPTR2X0=EPTR(300.-DELTA,100.,300.0)
c
c bottomside profile parameters .............................
C
      HMF1=HMF2
      HZ=HMF2
      HEF=HME
      B1=3.0
      B0 = B0_TAB(HOUR,SAX,SUX,NSEASN,RSSN,MODIP)
C!!!!!!! F1-REGION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      F1REG=.FALSE.
      HMF1=0.
      PNMF1=0.
      C1=0.  
      IF(.not.NIGHT) then
          FOF1=FOF1ED(ABSMBR,RSSN,XHI)
          IF(FOF1.gt.1.E-3) then
             F1REG=.TRUE.
             C1=.09+.11/DELA
             PNMF1=1.24E10*FOF1*FOF1
          endif
      endif
      NMF1=PNMF1
C!!!!!!! PARAMETER FOR E AND VALLEY-REGION !!!!!!!!!!!!!!!!!!!!!
      XDEL=XDELS(SEASON)/DELA
      DNDHBR=DNDS(SEASON)/DELA
      HDEEP=HPOL(HOUR,10.5/DELA,28.,SAX,SUX,1.,1.)
      WIDTH=HPOL(HOUR,17.8/DELA,45.+22./DELA,SAX,SUX,1.,1.)
      DEPTH=HPOL(HOUR,XDEL,81.,SAX,SUX,1.,1.)
      DLNDH=HPOL(HOUR,DNDHBR,.06,SAX,SUX,1.,1.)
      IF(DEPTH.ge.1.0) then
         IF(NIGHT) DEPTH=-DEPTH
         CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
         IF(EXT) then
            if(konsol.ge.0) WRITE(KONSOL,650)
	    width=0.
         endif
      endif
650   FORMAT(1X,'*NE* E-REGION VALLEY CAN NOT BE MODELLED')
      HEF=HME+WIDTH
      VNER = (1. - ABS(DEPTH) / 100.) * NME
c
c Parameters below E  .............................
c
      NMD=XMDED(XHI,RSSN,4.0E8)
      HMD=HPOL(HOUR,81.0,88.0,SAX,SUX,1.,1.)
      F(1)=HPOL(HOUR,0.02+0.03/DELA,0.05,SAX,SUX,1.,1.)
      F(2)=HPOL(HOUR,4.6,4.5,SAX,SUX,1.,1.)
      F(3)=HPOL(HOUR,-11.5,-4.0,SAX,SUX,1.,1.)
      FP1=F(1)
      FP2=-FP1*FP1/2.0
      FP30=(-F(2)*FP2-FP1+1.0/F(2))/(F(2)*F(2))
      FP3U=(-F(3)*FP2-FP1-1.0/F(3))/(F(3)*F(3))
c indermediate region between D and E region; parameters xkk
c and d1 are found such that the function reaches hdx/xdx/dxdh
      HDX=HMD+F(2)
      X=HDX-HMD
      XDX=NMD*EXP(X*(FP1+X*(FP2+X*FP30)))
      DXDX=XDX*(FP1+X*(2.0*FP2+X*3.0*FP30))
      X=HME-HDX
      XKK=-DXDX*X/(XDX*ALOG(XDX/NME))
c if exponent xkk is larger than xkkmax, then xkk will be set to xkkmax and
c d1 will be determined such that the point hdx/xdx is reached; derivative
c is no longer continuous.
      xkkmax=5.
      if(xkk.gt.xkkmax) then
         xkk=xkkmax
         d1=-alog(xdx/nme)/(x**xkk)
      else
         D1=DXDX/(XDX*XKK*X**(XKK-1.0))
      endif
C
C SEARCH FOR HMF1 ..................................................
C
924   IF(.not.F1REG) goto 380
         XE2H=XE2(HEF)
         CALL REGFA1(HEF,HMF2,XE2H,NMF2,0.001,NMF1,XE2,SCHALT,HMF1)
         IF(.not.SCHALT) GOTO 380
         if(konsol.ge.0) WRITE(KONSOL,11)
11       FORMAT(1X,'*NE* HMF1 IS NOT EVALUATED BY THE FUNCTION XE2')
         IREGFA=1
c
c change B1 and try again ..........................................
c
9244     IF(B1.GT.4.5) GOTO (7398,8922) IREGFA
         B1=B1+0.5
         if(konsol.ge.0) WRITE(KONSOL,902) B1-0.5,B1
902      FORMAT(6X,'CORR.: B1(OLD)=',F4.1,' B1(NEW)=',F4.1)
         GOTO 924
c
c omit F1 feature ....................................................
c
7398     if(konsol.ge.0) WRITE(KONSOL,9269)
9269     FORMAT(1X,'CORR.: NO F1 REGION, B1=3, C1=0.0')
         HMF1=0.
         NMF1=0.
         C1=0.0
         B1=3.
         F1REG=.FALSE.
380   continue
C
C SEARCH FOR HST [NE3(HST)=NME] ..........................................
C
      RRRR=0.5
      IF(F1REG) then
         hf1=hmf1
         xf1=nmf1
      else
         RATHH=0.5
3973     hf1=hef+(hmf2-hef)*RATHH
         xf1=xe3(hf1)
         IF(XF1.LT.NME) THEN
            RATHH=RATHH+.1
            GOTO 3973
         ENDIF
      endif
      h=hf1
      deh=10.
      XXMIN=XF1
      HHMIN=HF1
3895  h=h-deh
      if(h.lt.HEF) then
         h=h+2*deh
         deh=deh/10.
         if(deh.lt.1.) goto 3885
      endif
      XE3H=XE3(h)
      IF(XE3H.LT.XXMIN) then
         XXMIN=XE3H
         HHMIN=h
      endif
      if(XE3H.gt.NME) goto 3895
      CALL REGFA1(h,HF1,XE3H,XF1,0.001,NME,XE3,SCHALT,HST)
      STR=HST
      IF(.not.SCHALT) GOTO 360
3885  if(konsol.ge.0) WRITE(KONSOL,100)
100   FORMAT(1X,'*NE* HST IS NOT EVALUATED BY THE FUNCTION XE3')
      IREGFA=2
      IF(XXMIN/NME.LT.1.3) GOTO 9244
c
c assume linear interpolation between HZ and HEF ..................
c
8922  HZ=HHMIN+(HF1-HHMIN)*RRRR
      XNEHZ=XE3(HZ)
      if(xnehz-nme.lt.0.001) then
         RRRR=RRRR+.1
         GOTO 8922
      endif
      if(konsol.ge.0) WRITE(KONSOL,901) HZ,HEF
901   FORMAT(6X,'CORR.: LIN. APP. BETWEEN HZ=',F5.1,
     &          ' AND HEF=',F5.1)
      T=(XNEHZ-NME)/(HZ-HEF)
      HST=-333.
      GOTO 4933
c
c calculate HZ, D and T ............................................
c
360   HZ=(HST+HF1)/2.0
      D=HZ-HST
      T=D*D/(HZ-HEF-D)

C---------- CALCULATION OF NEUTRAL TEMPERATURE PARAMETER-------

4933  HTA=120.0
      HTE=3000.0
C
C CALCULATION FOR THE REQUIRED HEIGHT RANGE.......................
C
      IF(.NOT.F1REG) HMF1=HZ
      RETURN
        
* Error exits

8448  WRITE(*,8449) FILNAM
8449  FORMAT(' The file ',A30,'is not in your directory.')
      return
      END
