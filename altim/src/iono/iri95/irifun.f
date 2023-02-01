* "Light" version of the original IRIF13.F function in the
* IRI95 code.
C**************************************************************
C********** INTERNATIONAL REFERENCE IONOSPHERE ****************
C**************************************************************
C****************  FUNCTIONS,SUBROUTINES  *********************
C**************************************************************
C** IMPORTANT!! INITIALIZE (needs to be called before using
C**                 subroutines or functions)
C** NE:         XE1,ZERO,DXE1N,XE2,XE3,XE4,XE5,XE6,XE
C** PEAKS:      FOUT,XMOUT,HMF2ED,FOF1ED,FOEEDI,XMDED,GAMMA1
C** MAG. FIELD: GGM,FIELDG
C** FUNCTIONS:  REGFA1,TAL
C** TIME:       SOCO,HPOL,MODA
C** INTERPOL.:  B0POL,B0_TAB
C** EPSTEIN:    EPTR,EPST,EPSTEP,EPLA
C** INDICES:    TCON
C**************************************************************
C
C**************************************************************
C***  -------------------ADDRESSES------------------------  ***
C***  I  PROF. K. RAWER             DR. D. BILITZA       I  ***
C***  I  HERRENSTR. 43              GSFC CODE 933        I  ***
C***  I  7801 MARCH 1               GREENBELT MD 20771   I  ***
C***  I  F.R.G.                     USA                  I  ***
C***  ----------------------------------------------------  ***
C**************************************************************
C**************************************************************
C
C
C*************************************************************
C*************** ELECTRON DENSITY ****************************
C*************************************************************
C
C
      FUNCTION XE1(H)
c----------------------------------------------------------------
C REPRESENTING ELECTRON DENSITY(M-3) IN THE TOPSIDE IONOSPHERE
C (H=HMF2....1000 KM) BY HARMONIZED BENT-MODEL ADMITTING
C VARIABILITY OFGLOBAL PARAMETER ETA,ZETA,BETA,DELTA WITH
C GEOM. LATITUDE, SMOOTHED SOLAR FLUX AND CRITICAL FREQUENCY
C (SEE MAIN PROGRAM).
C [REF.:K.RAWER,S.RAMAKRISHNAN,1978]
c----------------------------------------------------------------
      COMMON    /BLOCK1/        HMF2,XNMF2,HMF1
      common    /BLO10/         BETA,ETA,DELTA,ZETA,EPTR1X0,EPTR2X0
      common    /ARGEXP/        ARGMAX
      dxdh = (1000.-hmf2)/700.
      x0 = 300. - delta
      xmx0 = (h-hmf2)/dxdh
      x = xmx0 + x0
      eptr1 = eptr(x,beta,394.5) - eptr1x0
      eptr2 = eptr(x,100.,300.0) - eptr2x0
      y = beta * eta * eptr1 + zeta * (100. * eptr2 - xmx0)
      y = y * dxdh
      if(abs(y).gt.argmax) y = sign(argmax,y)
      xe1 = xnmf2 * exp(-y)
      RETURN
      END
C
C
      REAL FUNCTION XE2(H)
C ELECTRON DENSITY FOR THE BOTTOMSIDE F-REGION (HMF1...HMF2).
      COMMON    /BLOCK1/HMF2,XNMF2,HMF1
      common    /BLOCK2/B0,B1,C1
      common	/ARGEXP/ARGMAX
      X=(HMF2-H)/B0
      if(x.le.0.0) x=0.0
      z=x**b1
      if(z.gt.argmax) z=argmax
      XE2=XNMF2*EXP(-z)/COSH(X)
      RETURN
      END
C
C
      REAL FUNCTION XE3(H)
C ELECTRON DENSITY FOR THE F1-LAYER (HZ.....HMF1).
      COMMON    /BLOCK1/        HMF2,XNMF2,HMF1
      common    /BLOCK2/        B0,B1,C1
      XE3=XE2(H)+XNMF2*C1*SQRT(ABS(HMF1-H)/B0)
      RETURN
      END
C
C
      REAL FUNCTION XE4(H)
C ELECTRON DENSITY FOR THE INDERMEDIUM REGION (HEF..HZ).
      COMMON    /BLOCK3/        HZ,T,HST,STR
      common    /BLOCK4/        HME,XNME,HEF
      IF(HST.GE.0.) THEN
         XE4=XE3(HZ+T/2.0-SIGN(1.0,T)*SQRT(T*(HZ-H+T/4.0)))
      ELSE
         XE4=XNME+T*(H-HEF)
      ENDIF
      RETURN
      END
C
C
      REAL FUNCTION XE5(H)
C ELECTRON DENSITY FOR THE E AND VALLEY REGION (HME..HEF).
      LOGICAL NIGHT
      COMMON    /BLOCK4/        HME,XNME,HEF
      common    /BLOCK5/        NIGHT,E(4)
      T3=H-HME
      T1=T3*T3*(E(1)+T3*(E(2)+T3*(E(3)+T3*E(4))))
      IF (NIGHT) THEN
         XE5=XNME*EXP(T1)
      ELSE
         XE5=XNME*(1+T1)
      ENDIF
      RETURN
      END
C
C
      REAL FUNCTION XE6(H)
C ELECTRON DENSITY FOR THE D REGION (HA...HME).
      COMMON    /BLOCK4/        HME,XNME,HEF
      common    /BLOCK6/        HMD,XNMD,HDX
      common    /BLOCK7/        D1,XKK,FP30,FP3U,FP1,FP2
      IF(H.LE.HDX) THEN
         Z=H-HMD
         FP3=FP3U
         IF(Z.GT.0.0) FP3=FP30
         XE6=XNMD*EXP(Z*(FP1+Z*(FP2+Z*FP3)))
      ELSE
         Z=HME-H
         XE6=XNME*EXP(-D1*Z**XKK)
      ENDIF
      RETURN
      END
C
C
      REAL FUNCTION XE(H)
C ELECTRON DENSITY BEETWEEN HA(KM) AND 1000 KM
C SUMMARIZING PROCEDURES  NE1....6;
      COMMON    /BLOCK1/HMF2,XNMF2,HMF1
      common    /BLOCK3/HZ,T,HST,STR
      common    /BLOCK4/HME,XNME,HEF
      IF (H.GE.HMF2) THEN
         XE=XE1(H)
      ELSE IF (H.GE.HMF1) THEN
         XE=XE2(H)
      ELSE IF (H.GE.HZ) THEN
         XE=XE3(H)
      ELSE IF (H.GE.HEF) THEN
         XE=XE4(H)
      ELSE IF (H.GE.HME) THEN
         XE=XE5(H)
      ELSE
         XE=XE6(H)
      ENDIF
      END
C
C*************************************************************
C************* PEAK VALUES ELECTRON DENSITY ******************
C*************************************************************
C
C
      real function FOUT(XMODIP,XLATI,XLONGI,UT,FF0)
C CALCULATES CRITICAL FREQUENCY FOF2/MHZ USING SUBROUTINE GAMMA1.
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME
C (DEC. HOURS), FF0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
      DIMENSION FF0(988)
      INTEGER QF(9)
      DATA QF/11,11,8,4,1,0,0,0,0/
      FOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,6,QF,9,76,13,988,FF0)
      RETURN
      END
C
C
      real function XMOUT(XMODIP,XLATI,XLONGI,UT,XM0)
C CALCULATES PROPAGATION FACTOR M3000 USING THE SUBROUTINE GAMMA1.
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME
C (DEC. HOURS), XM0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
      DIMENSION XM0(441)
      INTEGER QM(7)
      DATA QM/6,7,5,2,1,0,0/
      XMOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,4,QM,7,49,9,441,XM0)
      RETURN
      END
C
C
      REAL FUNCTION HMF2ED(XMAGBR,R,X,XM3)
C CALCULATES THE PEAK HEIGHT HMF2/KM FOR THE MAGNETIC
C LATITUDE XMAGBR/DEG. AND THE SMOOTHED ZUERICH SUNSPOT
C NUMBER R USING CCIR-M3000 XM3 AND THE RATIO X=FOF2/FOE.
C [REF. D.BILITZA ET AL., TELECOMM.J., 46, 549-553, 1979]
C D.BILITZA,1980.
      F1=(2.32E-3)*R+0.222
      F2=1.2-(1.16E-2)*EXP((2.39E-2)*R)
      F3=0.096*(R-25.0)/150.0
      DELM=F1*(1.0-R/150.0*EXP(-XMAGBR*XMAGBR/1600.0))/(X-F2)+F3
      HMF2ED=1490.0/(XM3+DELM)-176.0
      RETURN
      END
C
C
      REAL FUNCTION FOF1ED(YLATI,R,CHI)
c--------------------------------------------------------------
C CALCULATES THE F1 PEAK PLASMA FREQUENCY (FOF1/MHZ)
C FOR   DIP-LATITUDE (YLATI/DEGREE)
c       SMOOTHED ZURICH SUNSPOT NUMBER (R)
c       SOLAR ZENITH ANGLE (CHI/DEGREE)
C REFERENCE:
c       E.D.DUCHARME ET AL., RADIO SCIENCE 6, 369-378, 1971
C                                      AND 8, 837-839, 1973
c       HOWEVER WITH MAGNETIC DIP LATITUDE INSTEAD OF GEOMAGNETIC
c       DIPOLE LATITUDE, EYFRIG, 1979
C--------------------------------------------- D. BILITZA, 1988.
        COMMON/CONST/UMR
        FOF1ED = 0.0
        DLA =  YLATI
        CHI0 = 49.84733 + 0.349504 * DLA
        CHI100 = 38.96113 + 0.509932 * DLA
        CHIM = ( CHI0 + ( CHI100 - CHI0 ) * R / 100. )
        IF (CHI.GT.CHIM) RETURN
        F0 = 4.35 + DLA * ( 0.0058 - 1.2E-4 * DLA )
        F100 = 5.348 + DLA * ( 0.011 - 2.3E-4 * DLA )
        FS = F0 + ( F100 - F0 ) * R / 100.0
        XMUE = 0.093 + DLA * ( 0.0046 - 5.4E-5 * DLA ) + 3.0E-4 * R
        FOF1ED = FS * COS( CHI * UMR ) ** XMUE
        RETURN
        END
C
C
      REAL FUNCTION FOEEDI(COV,XHI,XHIM,XLATI)
C-------------------------------------------------------
C CALCULATES FOE/MHZ BY THE EDINBURGH-METHOD.
C INPUT: MEAN 10.7CM SOLAR RADIO FLUX (COV), GEOGRAPHIC
C LATITUDE (XLATI/DEG), SOLAR ZENITH ANGLE (XHI/DEG AND
C XHIM/DEG AT NOON).
C REFERENCE:
C       KOURIS-MUGGELETON, CCIR DOC. 6/3/07, 1973
C       TROST, J. GEOPHYS. RES. 84, 2736, 1979 (was used
C               to improve the nighttime varition)
C D.BILITZA--------------------------------- AUGUST 1986.
      COMMON/CONST/UMR
C variation with solar activity (factor A) ...............
      A=1.0+0.0094*(COV-66.0)
C variation with noon solar zenith angle (B) and with latitude (C)
      SL=COS(XLATI*UMR)
        IF(XLATI.LT.32.0) THEN
                SM=-1.93+1.92*SL
                C=23.0+116.0*SL
        ELSE
                SM=0.11-0.49*SL
                C=92.0+35.0*SL
        ENDIF
        if(XHIM.ge.90.) XHIM=89.999
        B = COS(XHIM*UMR) ** SM
C variation with solar zenith angle (D) ..........................
        IF(XLATI.GT.12.0) THEN
                SP=1.2
        ELSE
                SP=1.31
        ENDIF
C adjusted solar zenith angle during nighttime (XHIC) .............
      XHIC=XHI-3.*ALOG(1.+EXP((XHI-89.98)/3.))
      D=COS(XHIC*UMR)**SP
C determine foE**4 ................................................
      R4FOE=A*B*C*D
C minimum allowable foE (sqrt[SMIN])...............................
      SMIN=0.121+0.0015*(COV-60.)
      SMIN=SMIN*SMIN
      IF(R4FOE.LT.SMIN) R4FOE=SMIN
      FOEEDI=R4FOE**0.25
      RETURN
      END
C
C
      REAL FUNCTION XMDED(XHI,R,YW)
C D. BILITZA, 1978, CALCULATES ELECTRON DENSITY OF D MAXIMUM.
C XHI/DEG. IS SOLAR ZENITH ANGLE, R SMOOTHED ZURICH SUNSPOT NUMBER
C AND YW/M-3 THE ASSUMED CONSTANT NIGHT VALUE.
C [REF.: D.BILITZA, WORLD DATA CENTER A REPORT UAG-82,7,BOULDER,1981]
C corrected 4/25/97 - D. Bilitza
c
      COMMON/CONST/UMR
      if (xhi.lt.90) then
         Y = 6.05E8 + 0.088E8 * R
         yy = cos ( xhi * umr )
         ymd = y * exp( -0.1 / ( yy**2.7 ) )
         if (ymd.lt.yw) ymd = yw
	 xmded=ymd
      else
         XMDED=YW
      endif
      RETURN
      END
C
C
      REAL FUNCTION GAMMA1(SMODIP,SLAT,SLONG,HOUR,IHARM,NQ,
     &                          K1,M,MM,M3,SFE)
C CALCULATES GAMMA1=FOF2 OR M3000 USING CCIR NUMERICAL MAP
C COEFFICIENTS SFE(M3) FOR MODIFIED DIP LATITUDE (SMODIP/DEG)
C GEOGRAPHIC LATITUDE (SLAT/DEG) AND LONGITUDE (SLONG/DEG)
C AND UNIVERSIAL TIME (HOUR/DECIMAL HOURS).
C NQ(K1) IS AN INTEGER ARRAY GIVING THE HIGHEST DEGREES IN
C LATITUDE FOR EACH LONGITUDE HARMONIC.
C M=1+NQ1+2(NQ2+1)+2(NQ3+1)+... .
C SHEIKH,4.3.77.
      REAL*8 C(12),S(12),COEF(100),SUM
      DIMENSION NQ(K1),XSINX(13),SFE(M3)
      COMMON/CONST/UMR
      HOU=(15.0*HOUR-180.0)*UMR
      S(1)=SIN(HOU)
      C(1)=COS(HOU)
      DO I=2,IHARM
         C(I)=C(1)*C(I-1)-S(1)*S(I-1)
         S(I)=C(1)*S(I-1)+S(1)*C(I-1)
      ENDDO
      DO I=1,M
         MI=(I-1)*MM
         COEF(I)=SFE(MI+1)
         DO J=1,IHARM
            COEF(I)=COEF(I)+SFE(MI+2*J)*S(J)+SFE(MI+2*J+1)*C(J)
         ENDDO
      ENDDO
      SUM=COEF(1)
      SS=SIN(SMODIP*UMR)
      S3=SS
      XSINX(1)=1.0
      INDEX=NQ(1)
      DO J=1,INDEX
         SUM=SUM+COEF(1+J)*SS
         XSINX(J+1)=SS
         SS=SS*S3
      ENDDO
      XSINX(NQ(1)+2)=SS
      NP=NQ(1)+1
      SS=COS(SLAT*UMR)
      S3=SS
      DO J=2,K1
         S0=SLONG*(J-1.)*UMR
         S1=COS(S0)
         S2=SIN(S0)
         INDEX=NQ(J)+1
         DO L=1,INDEX
            NP=NP+1
            SUM=SUM+COEF(NP)*XSINX(L)*SS*S1
            NP=NP+1
            SUM=SUM+COEF(NP)*XSINX(L)*SS*S2
         ENDDO
         SS=SS*S3
      ENDDO
      GAMMA1=SUM
      RETURN
      END
C
C************************************************************
C*************** EARTH MAGNETIC FIELD ***********************
C**************************************************************
C
      SUBROUTINE GGM(ART,LONG,LATI,MLONG,MLAT)
C CALCULATES GEOMAGNETIC LONGITUDE (MLONG) AND LATITUDE (MLAT)
C FROM GEOGRAFIC LONGITUDE (LONG) AND LATITUDE (LATI) FOR ART=0
C AND REVERSE FOR ART=1. ALL ANGLES IN DEGREE.
C LATITUDE:-90 TO 90. LONGITUDE:0 TO 360 EAST.
      INTEGER ART
      REAL MLONG,MLAT,LONG,LATI
      COMMON/CONST/UMR
      ZPI=UMR*360.
      CBG=11.4*UMR
      CI=COS(CBG)
      SI=SIN(CBG)
      IF (ART.EQ.1) THEN
	 CBM=COS(MLAT*UMR)
	 SBM=SIN(MLAT*UMR)
	 CLM=COS(MLONG*UMR)
	 SLM=SIN(MLONG*UMR)
	 SBG=SBM*CI-CBM*CLM*SI
	 IF(ABS(SBG).GT.1.) SBG=SIGN(1.,SBG)
	 LATI=ASIN(SBG)
	 CBG=COS(LATI)
	 SLG=(CBM*SLM)/CBG
	 CLG=(SBM*SI+CBM*CLM*CI)/CBG
	 IF(ABS(CLG).GT.1.) CLG=SIGN(1.,CLG)
	 LONG=ACOS(CLG)
	 IF(SLG.LT.0.0) LONG=ZPI-LONG
	 LATI=LATI/UMR
	 LONG=LONG/UMR
	 LONG=LONG-69.8
	 IF(LONG.LT.0.0) LONG=LONG+360.0
      ELSE
	 YLG=LONG+69.8
	 CBG=COS(LATI*UMR)
	 SBG=SIN(LATI*UMR)
	 CLG=COS(YLG*UMR)
	 SLG=SIN(YLG*UMR)
	 SBM=SBG*CI+CBG*CLG*SI
	 IF(ABS(SBM).GT.1.) SBM=SIGN(1.,SBM)
	 MLAT=ASIN(SBM)
	 CBM=COS(MLAT)
	 SLM=(CBG*SLG)/CBM
	 CLM=(-SBG*SI+CBG*CLG*CI)/CBM
	 IF(ABS(CLM).GT.1.) CLM=SIGN(1.,CLM)
	 MLONG=ACOS(CLM)
	 IF(SLM.LT..0) MLONG=ZPI-MLONG
	 MLAT=MLAT/UMR
	 MLONG=MLONG/UMR
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE FIELDG(DLAT,DLONG,ALT,X,Y,Z,F,DIP,DEC,SMODIP)
C THIS IS A SPECIAL VERSION OF THE POGO 68/10 MAGNETIC FIELD
C LEGENDRE MODEL. TRANSFORMATION COEFF. G(144) VALID FOR 1973.
C INPUT: DLAT, DLONG=GEOGRAPHIC COORDINATES/DEG.(-90/90,0/360),
C        ALT=ALTITUDE/KM.
C OUTPUT: F TOTAL FIELD (GAUSS), Z DOWNWARD VERTICAL COMPONENT
C        X,Y COMPONENTS IN THE EQUATORIAL PLANE (X TO ZERO LONGITUDE).
C        DIP INCLINATION ANGLE(DEGREE). SMODIP RAWER'S MODFIED DIP.
C SHEIK,1977.
      DIMENSION H(144),XI(3),G(144),FEL1(72),FEL2(72)
      COMMON/CONST/UMR
      DATA FEL1/0.0, 0.1506723,0.0101742, -0.0286519, 0.0092606,
     & -0.0130846, 0.0089594, -0.0136808,-0.0001508, -0.0093977,
     & 0.0130650, 0.0020520, -0.0121956, -0.0023451, -0.0208555,
     & 0.0068416,-0.0142659, -0.0093322, -0.0021364, -0.0078910,
     & 0.0045586,  0.0128904, -0.0002951, -0.0237245,0.0289493,
     & 0.0074605, -0.0105741, -0.0005116, -0.0105732, -0.0058542,
     &0.0033268, 0.0078164,0.0211234, 0.0099309, 0.0362792,
     &-0.0201070,-0.0046350,-0.0058722,0.0011147,-0.0013949,
     & -0.0108838,  0.0322263, -0.0147390,  0.0031247, 0.0111986,
     & -0.0109394,0.0058112,  0.2739046, -0.0155682, -0.0253272,
     &  0.0163782, 0.0205730,  0.0022081, 0.0112749,-0.0098427,
     & 0.0072705, 0.0195189, -0.0081132, -0.0071889, -0.0579970,
     & -0.0856642, 0.1884260,-0.7391512, 0.1210288, -0.0241888,
     & -0.0052464, -0.0096312, -0.0044834, 0.0201764,  0.0258343,
     &0.0083033,  0.0077187/
      DATA FEL2/0.0586055,0.0102236,-0.0396107,
     & -0.0167860, -0.2019911, -0.5810815,0.0379916,  3.7508268,
     & 1.8133030, -0.0564250, -0.0557352, 0.1335347, -0.0142641,
     & -0.1024618,0.0970994, -0.0751830,-0.1274948, 0.0402073,
     &  0.0386290, 0.1883088,  0.1838960, -0.7848989,0.7591817,
     & -0.9302389,-0.8560960, 0.6633250, -4.6363869, -13.2599277,
     & 0.1002136,  0.0855714,-0.0991981, -0.0765378,-0.0455264,
     &  0.1169326, -0.2604067, 0.1800076, -0.2223685, -0.6347679,
     &0.5334222, -0.3459502,-0.1573697,  0.8589464, 1.7815990,
     &-6.3347645, -3.1513653, -9.9927750,13.3327637, -35.4897308,
     &37.3466339, -0.5257398,  0.0571474, -0.5421217,  0.2404770,
     & -0.1747774,-0.3433644, 0.4829708,0.3935944, 0.4885033,
     &  0.8488121, -0.7640999, -1.8884945, 3.2930784,-7.3497229,
     & 0.1672821,-0.2306652, 10.5782146, 12.6031065, 8.6579742,
     & 215.5209961, -27.1419220,22.3405762,1108.6394043/
      K=0
      DO 10 I=1,72
      K=K+1
      G(K)=FEL1(I)
10    G(72+K)=FEL2(I)
      RLAT=DLAT*UMR
      CT=SIN(RLAT)
      ST=COS(RLAT)
      NMAX=11
      D=SQRT(40680925.0-272336.0*CT*CT)
      RLONG=DLONG*UMR
      CP=COS(RLONG)
      SP=SIN(RLONG)
      ZZZ=(ALT+40408589.0/D)*CT/6371.2
      RHO=(ALT+40680925.0/D)*ST/6371.2
      XXX=RHO*CP
      YYY=RHO*SP
      RQ=1.0/(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      XI(1)=XXX*RQ
      XI(2)=YYY*RQ
      XI(3)=ZZZ*RQ
      IHMAX=NMAX*NMAX+1
      LAST=IHMAX+NMAX+NMAX
      IMAX=NMAX+NMAX-1
      DO I=IHMAX,LAST
         H(I)=G(I)
      ENDDO
      DO K=1,3,2
	 I=IMAX
	 IH=IHMAX
300      IL=IH-I
	 F1=2./(I-K+2.)
	 X1=XI(1)*F1
	 Y1=XI(2)*F1
	 Z1=XI(3)*(F1+F1)
	 I=I-2
	 IF((I-1).LT.0) GOTO 400
	 IF((I-1).EQ.0) GOTO 500
	 DO M=3,I,2
	    H(IL+M+1)=G(IL+M+1)+Z1*H(IH+M+1)+X1*(H(IH+M+3)-H(IH+M-1))-
     &	    Y1*(H(IH+M+2)+H(IH+M-2))
	    H(IL+M)=G(IL+M)+Z1*H(IH+M)+X1*(H(IH+M+2)-H(IH+M-2))+
     &	    Y1*(H(IH+M+3)+H(IH+M-1))
         ENDDO
500      H(IL+2)=G(IL+2)+Z1*H(IH+2)+X1*H(IH+4)-Y1*(H(IH+3)+H(IH))
         H(IL+1)=G(IL+1)+Z1*H(IH+1)+Y1*H(IH+4)+X1*(H(IH+3)-H(IH))
400      H(IL)=G(IL)+Z1*H(IH)+2.0*(X1*H(IH+1)+Y1*H(IH+2))
         IH=IL
         IF(I.GE.K) GOTO 300
      ENDDO
      S=0.5*H(1)+2.0*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      XT=(RQ+RQ)*SQRT(RQ)
      X=XT*(H(3)-S*XXX)
      Y=XT*(H(4)-S*YYY)
      Z=XT*(H(2)-S*ZZZ)
      F=SQRT(X*X+Y*Y+Z*Z)
      BRH0=Y*SP+X*CP
      Y=Y*CP-X*SP
      X=Z*ST-BRH0*CT
      Z=-Z*CT-BRH0*ST
      zdivf=z/f
      IF(ABS(zdivf).GT.1.) zdivf=SIGN(1.,zdivf)
      DIP=ASIN(zdivf)
      ydivs=y/sqrt(x*x+y*y)
      IF(ABS(ydivs).GT.1.) ydivs=SIGN(1.,ydivs)
      DEC=ASIN(ydivs)
      dipdiv=DIP/SQRT(DIP*DIP+ST)
      IF(ABS(dipdiv).GT.1.) dipdiv=SIGN(1.,dipdiv)
      SMODIP=ASIN(dipdiv)
      DIP=DIP/UMR
      DEC=DEC/UMR
      SMODIP=SMODIP/UMR
      RETURN
      END
C
C
C************************************************************
C*********** INTERPOLATION AND REST ***************************
C**************************************************************
C
C
      SUBROUTINE REGFA1(X11,X22,FX11,FX22,EPS,FW,F,SCHALT,X)
C REGULA-FALSI-PROCEDURE TO FIND X WITH F(X)-FW=0. X1,X2 ARE THE
C STARTING VALUES. THE COMUTATION ENDS WHEN THE X-INTERVAL
C HAS BECOME LESS THAN EPS . IF SIGN(F(X1)-FW)= SIGN(F(X2)-FW)
C THEN SCHALT=.TRUE.
      LOGICAL L1,LINKS,K,SCHALT
      SCHALT=.FALSE.
      EP=EPS
      X1=X11
      X2=X22
      F1=FX11-FW
      F2=FX22-FW
      K=.FALSE.
      L1=.FALSE.
      NG=2
      LFD=0
      IF(F1*F2.LE.0.0) GOTO 200
        X=0.0
        SCHALT=.TRUE.
        RETURN
200   X=(X1*F2-X2*F1)/(F2-F1)
      GOTO 400
300     L1=LINKS
        DX=(X2-X1)/NG
        IF(.NOT.LINKS) DX=DX*(NG-1)
        X=X1+DX
400   FX=F(X)-FW
      LFD=LFD+1
      IF(LFD.GT.20) THEN
        EP=EP*10.
        LFD=0
      ENDIF
      LINKS=(F1*FX.GT.0.0)
      K=.NOT.K
      IF(LINKS) THEN
        X1=X
        F1=FX
      ELSE
        X2=X
        F2=FX
      ENDIF
      IF(ABS(X2-X1).LE.EP) GOTO 800
      IF(K) GOTO 300
      IF(LINKS.NEQV.L1) NG=2*NG
      GOTO 200
800   RETURN
      END
C
C
      SUBROUTINE TAL(SHABR,SDELTA,SHBR,SDTDH0,AUS6,SPT)
C CALCULATES THE COEFFICIENTS SPT FOR THE POLYNOMIAL
C Y(X)=1+SPT(1)*X**2+SPT(2)*X**3+SPT(3)*X**4+SPT(4)*X**5
C TO FIT THE VALLEY IN Y, REPRESENTED BY:
C Y(X=0)=1, THE X VALUE OF THE DEEPEST VALLEY POINT (SHABR),
C THE PRECENTAGE DEPTH (SDELTA), THE WIDTH (SHBR) AND THE
C DERIVATIVE DY/DX AT THE UPPER VALLEY BOUNDRY (SDTDH0).
C IF THERE IS AN UNWANTED ADDITIONAL EXTREMUM IN THE VALLEY
C REGION, THEN AUS6=.TRUE., ELSE AUS6=.FALSE..
C FOR -SDELTA THE COEFF. ARE CALCULATED FOR THE FUNCTION
C Y(X)=EXP(SPT(1)*X**2+...+SPT(4)*X**5).
      DIMENSION SPT(4)
      LOGICAL AUS6
      Z1=-SDELTA/(100.0*SHABR*SHABR)
      IF(SDELTA.GT.0.) GOTO 500
      SDELTA=-SDELTA
      Z1=ALOG(1.-SDELTA/100.)/(SHABR*SHABR)
500   Z3=SDTDH0/(2.*SHBR)
      Z4=SHABR-SHBR
      SPT(4)=2.0*(Z1*(SHBR-2.0*SHABR)*SHBR+Z3*Z4*SHABR)/
     &  (SHABR*SHBR*Z4*Z4*Z4)
      SPT(3)=Z1*(2.0*SHBR-3.0*SHABR)/(SHABR*Z4*Z4)-
     &  (2.*SHABR+SHBR)*SPT(4)
      SPT(2)=-2.0*Z1/SHABR-2.0*SHABR*SPT(3)-3.0*SHABR*SHABR*SPT(4)
      SPT(1)=Z1-SHABR*(SPT(2)+SHABR*(SPT(3)+SHABR*SPT(4)))
      AUS6=.FALSE.
      B=4.*SPT(3)/(5.*SPT(4))+SHABR
      C=-2.*SPT(1)/(5*SPT(4)*SHABR)
      Z2=B*B/4.-C
      IF(Z2.LT.0.0) GOTO 300
      Z3=SQRT(Z2)
      Z1=B/2.
      Z2=-Z1+Z3
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.
      IF (ABS(Z3).GT.1.E-15) GOTO 400
      Z2=C/Z2
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.
      RETURN
400   Z2=-Z1-Z3
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.
300   RETURN
      END
C
C
C******************************************************************
C********** ZENITH ANGLE, DAY OF YEAR, TIME ***********************
C******************************************************************
C
C
      subroutine soco (ld,t,flat,Elon,
     &        DECLIN, ZENITH, SUNRSE, SUNSET)
c--------------------------------------------------------------------
c       s/r to calculate the solar declination, zenith angle, and
c       sunrise & sunset times  - based on Newbern Smith's algorithm
c       [leo mcnamara, 1-sep-86, last modified 16-jun-87]
c       {dieter bilitza, 30-oct-89, modified for IRI application}
c
c in:   ld      local day of year
c       t       local hour (decimal)
c       flat    northern latitude in degrees
c       elon    east longitude in degrees
c
c out:  declin      declination of the sun in degrees
c       zenith      zenith angle of the sun in degrees
c       sunrse      local time of sunrise in hours
c       sunset      local time of sunset in hours
c-------------------------------------------------------------------
c
      common /const/umr
      common /const1/humr
c ampludes of Fourier coefficients  --  1955 epoch.................
      data    p1,p2,p3,p4,p6 /
     &0.017203534,0.034407068,0.051610602,0.068814136,0.103221204 /
c
c s/r is formulated in terms of WEST longitude.......................
      wlon = 360. - Elon
c
c time of equinox for 1980...........................................
      td = ld + (t + Wlon/15.) / 24.
      te = td + 0.9369
c
c declination of the sun..............................................
      dcl = 23.256 * sin(p1*(te-82.242)) + 0.381 * sin(p2*(te-44.855))
     &    + 0.167 * sin(p3*(te-23.355)) - 0.013 * sin(p4*(te+11.97))
     &    + 0.011 * sin(p6*(te-10.41)) + 0.339137
      DECLIN = dcl
      dc = dcl * umr
c
c the equation of time................................................
      tf = te - 0.5
      eqt = -7.38*sin(p1*(tf-4.)) - 9.87*sin(p2*(tf+9.))
     &    + 0.27*sin(p3*(tf-53.)) - 0.2*cos(p4*(tf-17.))
      et = eqt * umr / 4.
c
      fa = flat * umr
      phi = humr * ( t - 12.) + et
c
      a = sin(fa) * sin(dc)
      b = cos(fa) * cos(dc)
      cosx = a + b * cos(phi)
      if(abs(cosx).gt.1.) cosx=sign(1.,cosx)
      zenith = acos(cosx) / umr
c
c calculate sunrise and sunset times --  at the ground...........
c see Explanatory Supplement to the Ephemeris (1961) pg 401......
c sunrise at height h metres is at...............................
c       chi(h) = 90.83 + 0.0347 * sqrt(h)........................
c this includes corrections for horizontal refraction and........
c semi-diameter of the solar disk................................
      ch = cos(90.83 * umr)
      cosphi = (ch -a ) / b
c if abs(secphi) > 1., sun does not rise/set.....................
c allow for sun never setting - high latitude summer.............
      secphi = 999999.
      if(cosphi.ne.0.) secphi = 1./cosphi
      sunset = 99.
      sunrse = 99.
      if(secphi.gt.-1.0.and.secphi.le.0.) return
c allow for sun never rising - high latitude winter..............
      sunset = -99.
      sunrse = -99.
      if(secphi.gt.0.0.and.secphi.lt.1.) return
c
      if(cosphi.gt.1.) cosphi=sign(1.,cosphi)
      phi = acos(cosphi)
      et = et / humr
      phi = phi / humr
      sunrse = 12. - phi - et
      sunset = 12. + phi - et
      if(sunrse.lt.0.) sunrse = sunrse + 24.
      if(sunset.ge.24.) sunset = sunset - 24.
c
      return
      end
c
C
      FUNCTION HPOL(HOUR,TW,XNW,SA,SU,DSA,DSU)
C-------------------------------------------------------
C PROCEDURE FOR SMOOTH TIME-INTERPOLATION USING EPSTEIN
C STEP FUNCTION AT SUNRISE (SA) AND SUNSET (SU). THE
C STEP-WIDTH FOR SUNRISE IS DSA AND FOR SUNSET DSU.
C TW,NW ARE THE DAY AND NIGHT VALUE OF THE PARAMETER TO
C BE INTERPOLATED. SA AND SU ARE TIME OF SUNRIES AND
C SUNSET IN DECIMAL HOURS.
C BILITZA----------------------------------------- 1979.
      IF (SU.LT.-25.) THEN
         HPOL=XNW
      ELSE IF (SU.GT.25.) THEN
         HPOL=TW
      ELSE
         HPOL=XNW+(TW-XNW)*EPST(HOUR,DSA,SA)+
     &  (XNW-TW)*EPST(HOUR,DSU,SU)
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE MODA(IN,IYEAR,MONTH,IDAY,IDOY,NRDAYMO)
C-------------------------------------------------------------------
C CALCULATES DAY OF YEAR (IDOY, ddd) FROM YEAR (IYEAR, yy or yyyy),
C MONTH (MONTH, mm) AND DAY OF MONTH (IDAY, dd) IF IN=0, OR MONTH
C AND DAY FROM YEAR AND DAY OF YEAR IF IN=1. NRDAYMO is an output
C parameter providing the number of days in the specific month.
C-------------------------------------------------------------------
      DIMENSION       MM(12)
      DATA            MM/31,28,31,30,31,30,31,31,30,31,30,31/

      IMO=0
      MOBE=0
c
c  leap year rule: years evenly divisible by 4 are leap years, except
c  years also evenly divisible by 100 are not leap years, except years also
c  evenly divisible by 400 are leap years. The year 2000 therefore is a
c  leap year. The 100 and 400 year exception rule
c       if((iyear/4*4.eq.iyear).and.(iyear/100*100.ne.iyear)) mm(2)=29
c  will become important again in the year 2100 which is not a leap year.
c
      mm(2)=28
      if(iyear/4*4.eq.iyear) mm(2)=29

      IF(IN.le.0) then
         mosum=0
         do i=1,month-1
            mosum=mosum+mm(i)
         enddo
         idoy=mosum+iday
         nrdaymo=mm(month)
      else
5        IMO=IMO+1
         IF(IMO.GT.12) return
         MOOLD=MOBE
         nrdaymo=mm(imo)
         MOBE=MOBE+nrdaymo
         IF(MOBE.LT.IDOY) GOTO 5
         MONTH=IMO
         IDAY=IDOY-MOOLD
      endif
      RETURN
      END
c
c
      REAL FUNCTION B0_TAB ( HOUR, SAX, SUX, NSEASN, R, ZMODIP)
C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON is northern season with
C ISEASON=1 northern spring), low and high solar activity Rz12=10,
C 100 (IR=1,2), and low and middle modified dip latitudes 18 and 45
C degress (ILATI=1,2). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
C corrected to include a smooth transition at the modip equator
C and no discontinuity at the equatorial change in season.
C JAN 1993 ---------------------------------------- Dieter Bilitza
C
      REAL      NITVAL
      DIMENSION B0F(2,4,2,2),bfr(2,2,2),bfd(2,2),zx(5),g(6),dd(5)
      DATA      B0F/114.,64.0,134.,77.0,128.,66.0,75.,73.0,
     &              113.,115.,150.,116.,138.,123.,94.,132.,
     &              72.0,84.0,83.0,89.0,75.0,85.0,57.,76.0,
     &              102.,100.,120.,110.,107.,103.,76.,86.0/
      data    zx/45.,70.,90.,110.,135./,dd/2.5,2.,2.,2.,2.5/

C jseasn is southern hemisphere season
      jseasn=nseasn+2
      if(jseasn.gt.4) jseasn=jseasn-4

      zz = zmodip + 90.
      zz0 = 0.

C Interpolation in Rz12: linear from 10 to 100
      DO ISL=1,2
         DO ISD=1,2
            bfr(isd,1,isl) = b0f(isd,nseasn,1,isl) +
     &      (b0f(isd,nseasn,2,isl) - b0f(isd,nseasn,1,isl))/90.*(R-10.)
            bfr(isd,2,isl) = b0f(isd,jseasn,1,isl) +
     &      (b0f(isd,jseasn,2,isl) - b0f(isd,jseasn,1,isl))/90.*(R-10.)
         enddo

C Interpolation day/night with transitions at SAX (sunrise) and SUX (sunset)
         do iss=1,2
            DAYVAL = BFR(1,ISS,ISL)
            NITVAL = BFR(2,ISS,ISL)
            BFD(iss,ISL) = HPOL(HOUR,DAYVAL,NITVAL,SAX,SUX,1.,1.)
         enddo
      enddo

C Interpolation with epstein-transitions in modified dip latitude.
C Transitions at +/-18 and +/-45 degrees; constant above +/-45.
C
C g(1:5) are the latitudinal slopes; g(1) is for the region from -90
C to -45 degrees, g(2) for -45/-20, g(3) for -20/0, g(4) for 0/20,
C g(5) for 20/45, and g(6) for 45/90. B0=bfd(2,2) at modip = -90,
C bfd(2,2) at modip = -45, bfd(2,1) at modip = -20, bfd(2,1)+delta at
C modip = -10 and 0, bfd(1,1) at modip = 20, bfd(1,2) at modip = 45 and 90.

      g(1) = 0.
      g(2) = ( bfd(2,1) - bfd(2,2) ) / 25.
      g(5) = ( bfd(1,2) - bfd(1,1) ) / 25.
      g(6) = 0.
      if(bfd(2,1).gt.bfd(1,1)) then
         g(3) = g(2) / 4.
         yb4 = bfd(2,1) + 20. * g(3)
         g(4) = ( bfd(1,1) - yb4 ) / 20.
      else
         g(4) = g(5) / 4.
         yb5 = bfd(1,1) - 20. * g(4)
         g(3) = ( yb5 - bfd(2,1) ) / 20.
      endif
      bb0 = bfd(2,2)
      SUM = bb0
      DO I=1,5
         aa = eptr(zz ,dd(i),zx(i))
         bb = eptr(zz0,dd(i),zx(i))
         DSUM = (G(I+1) - G(I)) * (AA-BB) * dd(i)
         SUM = SUM + DSUM
      enddo
      B0_TAB = SUM
      RETURN
      END
c
C
C *********************************************************************
C ************************ EPSTEIN FUNCTIONS **************************
C *********************************************************************
C REF:  H. G. BOOKER, J. ATMOS. TERR. PHYS. 39, 619-623, 1977
C       K. RAWER, ADV. SPACE RES. 4, #1, 11-15, 1984
C *********************************************************************
C
C
      REAL FUNCTION EPTR ( X, SC, HX )
C ------------------------------------------------------------ TRANSITION
      COMMON/ARGEXP/ARGMAX
      D1 = ( X - HX ) / SC
      IF (D1.LE.-ARGMAX) THEN
         EPTR = 0
      ELSE IF (D1.GT.ARGMAX) THEN
         EPTR = D1
      ELSE
         EPTR = ALOG ( 1. + EXP( D1 ))
      ENDIF
      RETURN
      END
C
C
      REAL FUNCTION EPST ( X, SC, HX )
C -------------------------------------------------------------- STEP
      COMMON/ARGEXP/ARGMAX
      D1 = ( X - HX ) / SC
      IF (D1.LE.-ARGMAX) THEN
         EPST = 0.
      ELSE IF (D1.GT.ARGMAX) THEN
         EPST = 1.
      ELSE
         EPST = 1. / ( 1. + EXP( -D1 ))
      ENDIF
      RETURN
      END
C
C
      subroutine tcon(yr,mm,day,idn,rz,ig,rsn,nmonth)
c----------------------------------------------------------------
c input:        yr,mm,day       year(yyyy),month(mm),day(dd)
c               idn             day of year(ddd)
c output:       rz(3)           12-month-smoothed solar sunspot number
c               ig(3)           12-month-smoothed IG index
c               rsn             interpolation parameter
c               nmonth          previous or following month depending
c                               on day
c
c rz(1), ig(1) contain the indices for the month mm and rz(2), ig(2)
c for the previous month (if day less than 15) or for the following
c month (otherwise). These indices are for the mid of the month. The
c indices for the given day are obtained by linear interpolation and
c are stored in rz(3) and ig(3).
c
      character*80 filenm
      integer      yr, mm, day, iflag, iyst, iyend,iymst
      integer      imst,iymend,unit,freeunit
      real         ionoindx(1202),indrz(1202)
      real         ig(3),rz(3)
      save         ionoindx,indrz,iflag,iyst,iymst,iymend,imst
c
c Rz12 and IG are determined from the file IG_RZ.DAT which has the
c following structure:
c day, month, year of the last update of this file,
c start month, start year, end month, end year,
c the 12 IG indices (13-months running mean) for the first year,
c the 12 IG indices for the second year and so on until the end year,
c the 12 Rz indices (13-months running mean) for the first year,
c the 12 Rz indices for the second year and so on until the end year.
c The inteporlation procedure also requires the IG and Rz values for
c the month preceeding the start month and the IG and Rz values for the
c month following the end month. These values are also included in IG_RZ.
c
c A negative Rz index means that the given index is the 13-months-running
c mean of the solar radio flux (F10.7). The close correlation between (Rz)12
c and (F10.7)12 is used to derive the (Rz)12 indices.
c
c An IG index of -111 indicates that no IG values are available for the
c time period. In this case a correlation function between (IG)12 and (Rz)12
c is used to obtain (IG)12.
c
c The computation of the 13-month-running mean for month M requires the indices
c for the six months preceeding M and the six months following M (month: M-6,
c ..., M+6). To calculate the current running mean one therefore requires
c predictions of the indix for the next six months. Starting from six months
c before the UPDATE DATE (listed at the top of the file) and onward the
c indices are therefore based on indices predictions.
c
      if(iflag.eq.0) then
	 unit=freeunit()
	 filenm='/user/altim'
	 call checkenv('ALTIM',filenm,l)
	 filenm(l+1:)="/data/indices/ig_rz.dat"
         open(unit=unit,file=filenm,status='old')
c Read the update date, the start date and the end date (mm,yyyy), and
c get number of data points to read.
         read(unit,*) iupd,iupm,iupy
         read(unit,*) imst,iyst, imend, iyend
         iymst=iyst*100+imst
         iymend=iyend*100+imend
c inum_vals= 12-imst+1+(iyend-iyst-1)*12 +imend + 2
c            1st year \ full years       \last y\ before & after
         inum_vals= 3-imst+(iyend-iyst)*12 +imend
c Read all the ionoindx and indrz values
         read(unit,*) (ionoindx(i),i=1,inum_vals)
         read(unit,*) (indrz(i),i=1,inum_vals)
         do jj=1,inum_vals
            rrr=indrz(jj)
            if(rrr.lt.0.0) then
               covr=abs(rrr)
               rrr=33.52*sqrt(covr+85.12)-408.99
               if(rrr.lt.0.0) rrr=0.0
               indrz(jj)=rrr
            endif
            if(ionoindx(jj).le.-90.) then
               zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
               if(zi.gt.274.0) zi=274.0
               ionoindx(jj)=zi
	    endif
         enddo
         close(unit)
         iflag = 1
      endif

      iytmp=yr*100+mm
      if (iytmp .lt. iymst .or. iytmp .gt. iymend) then
         write(6,8000) iytmp,iymst,iymend
 8000    format(1x,I10,'** OUT OF RANGE **'/,5x,
     &  'The file IG_RZ.DAT which contains the indices Rz12',
     &  ' and IG12'/5x,'currently only covers the time period',
     &  ' (yymm) : ',I6,'-',I6)
         nmonth=-1
         return
      endif

c     num=12-imst+1+(yr-iyst-1)*12+mm+1
      num=2-imst+(yr-iyst)*12+mm

      rz(1)=indrz(num)
      ig(1)=ionoindx(num)
      midm=15
      if(mm.eq.2) midm=14
      call MODA(0,yr,mm,midm,idd1,nrdaym)
      if(day.ge.midm) then
         imm2=mm+1
         if(imm2.gt.12) then
            imm2=1
            iyy2=yr+1
            idd2=380
            if(yr/4*4.eq.yr) idd2=381
         else
            iyy2=yr
            midm=15
            if(imm2.eq.2) midm=14
            call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
         endif
         rz(2)=indrz(num+1)
         ig(2)=ionoindx(num+1)
         rsn=(idn-idd1)*1./(idd2-idd1)
         rz(3)=rz(1)+(rz(2)-rz(1))*rsn
         ig(3)=ig(1)+(ig(2)-ig(1))*rsn
      else
         imm2=mm-1
         if(imm2.lt.1) then
            imm2=12
            idd2=-16
            iyy2=yr-1
         else
            iyy2=yr
            midm=15
            if(imm2.eq.2) midm=14
            call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
         endif
         rz(2)=indrz(num-1)
         ig(2)=ionoindx(num-1)
         rsn=(idn-idd2)*1./(idd1-idd2)
         rz(3)=rz(2)+(rz(1)-rz(2))*rsn
         ig(3)=ig(2)+(ig(1)-ig(2))*rsn
      endif
      nmonth=imm2
      return
      end
