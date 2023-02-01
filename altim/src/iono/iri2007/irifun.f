c irifun.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
C
C Functions and subroutines for the International Reference
C Ionosphere model. These functions and subroutines are called by
C IRI_SUB subroutine (IRISUB.FOR).
c
c-----------------------------------------------------------------------
C
c i/o units:  6   messages (during execution) to monitor
c
c-----------------------------------------------------------------------
c changes from IRIFU9 to IRIF10:
c       SOCO for solar zenith angle
c       ACOS and ASIN argument forced to be within -1 / +1
c       EPSTEIN functions corrected for large arguments
c-----------------------------------------------------------------------
c changes from IRIF10 to IRIF11:
c       LAY subroutines introduced
c       TEBA corrected for 1400 km
c-----------------------------------------------------------------------
c changes from IRIF11 to IRIF12:
C       Neutral temperature subroutines now in CIRA86.FOR
C       TEDER changed
C       All names with 6 or more characters replaced
C       10/29/91 XEN: 10^ in loop, instead of at the end
C       1/21/93 B0_TAB instead of B0POL
C       9/22/94 Alleviate underflow condition in IONCOM exp()
c-----------------------------------------------------------------------
c changes from IRIF12 to IRIF13:
C        9/18/95 MODA: add leap year and number of days in month
C        9/29/95 replace F2out with FOUT and XMOUT.
C       10/ 5/95 add TN and DTNDH; earlier in CIRA86.FOR
C       10/ 6/95 add TCON for reading indices
C       10/20/95 MODA: IN=1 MONTH=IMO
C       10/20/95 TCON: now includes RZ interpolation
C       11/05/95 IONCOM->IONCO1, added IONCOM_new, IONCO2
C       11/05/95 LSTID added for strom-time updating
C       11/06/95 ROGUL: transition 20. instead of 15.
C       12/01/95 add UT_LT for (date-)correct UT<->LT conversion
C       01/16/96 TCON: add IMST to SAVE statement
C       02/02/96 ROGUL: 15. reinstated
C       02/07/96 UT_LT: ddd, dddend integer, no leap year 2000
C       03/15/96 ZERO: finding delta for topside
C       03/18/96 UT_LT: mode=1, change of year
C       12/09/96 since 2000 is leap, delete y/100*100 condition
C       04/25/97 XMDED: minimal value also daytime
C       05/18/98 TCON: changes to IG_RZ (update date); -R = Cov
C       05/19/98 Replaced IONCO2&APROK; HEI,XHI in IONCOM_NEW
C       10/01/98 added INITIALIZE
C       04/30/99 MODA: reset bb(2)=28
C       11/08/99 avoid negative x value in function XE2. Set x=0.
C       11/30/99 EXIT in APROK replaced with GOTO (N. Smirnova)
c-----------------------------------------------------------------------
c changes from IRIF13 to IRIFUN:
C-Version-mm/dd/yy-description [person reporting correction]
C 2000.01 05/09/00 Include B0_98 subroutine to replace B0_NEW and B0POL
C 2000.02 05/18/00 Include Elteik and spharm_ik for Te
C 2000.03 06/09/00 Include XE3, XE4, XE
C 2000.04 06/11/00 Include f1_c1, f1_prob, modified fof1ed
C 2000.05 10/30/00 Include vdrift
C 2000.06 04/15/01 Include IGRF_SUB subroutine for IK Te model
C 2000.07 05/07/01 Include storm subroutine STORM and Ap access s/w
C 2000.08 09/07/01 APF: if(j1.eq.j2) -> if(IY.eq.j2) [P. Wilkinson]
C 2000.09 09/07/01 CONVER: LO2 = MOD(LO1,20)+1 [P. Webb,D. Pesnell]
C 2000.10 02/20/02 CONVER/DATA: 105.78 -> 015.78 [A. Shovkoplyas]
C 2000.11 10/28/02 replace TAB/6 blanks, enforce 72/line [D. Simpson]
C 2000.12 11/08/02 removing unused variables (corr); apf0 removed
C 2000.13 11/26/02 apf() using keyed access to ap.dat file; apf->apf1
C 2000.14 11/27/02 changed F1_PROB; always 6 preceeding spaces
C 2005.01 03/09/05 CALION,INVDPC,CALNE for new Ne, Ni models
C 2005.01 11/14/05 APF_ONLY for F107D;
C 2005.01 11/14/05 spreadf_brazil; added constraint 0<=P<=1
C 2005.02 05/11/06 NeQuick: XE1,TOPQ, M3000HM; stormvd,
C 2005.02 03/27/07 STORM: hourly interpolation of Ap  [A. Oinats]
C 2007.00 05/18/07 Release of IRI-2007
C 2007.01 09/19/07 vdrift et al.: without *8 (no change in results)
C
C**************************************************************
C********** INTERNATIONAL REFERENCE IONOSPHERE ****************
C**************************************************************
C****************  FUNCTIONS,SUBROUTINES  *********************
C**************************************************************
C** NE:         XE1,XE2,XE3,XE4,XE5,XE6,XE
C** PEAKS:      FOUT,XMOUT,HMF2ED,FOF1ED,f1_c1,f1_prob,FOEEDI,XMDED,
C**		          GAMMA1
C** PROFILE PAR:B0_98,TAL
C** MAG. FIELD: GGM,FIELDG,CONVER(Geom. Corrected Latitude)
C** FUNCTIONS:  REGFA1
C** TIME:       SOCO,HPOL,MODA,UT_LT
C** EPSTEIN:    EPTR,EPST,EPSTEP,EPLA
C** INDICES:    TCON,APF,APF_ONLY
C** Storm:   	LSTID,storm
C**************************************************************
C
C**************************************************************
C***  -------------------ADDRESSES------------------------  ***
C***  I  PROF. K. RAWER             DR. D. BILITZA       I  ***
C***  I  HERRENSTR. 43              GSFC CODE 632        I  ***
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

        FUNCTION XE1(H)
c----------------------------------------------------------------
C REPRESENTING ELECTRON DENSITY(M-3) IN THE TOPSIDE IONOSPHERE
C (H=HMF2....1000 KM) BY HARMONIZED BENT-MODEL ADMITTING
C VARIABILITY OFGLOBAL PARAMETER ETA,ZETA,BETA,DELTA WITH
C GEOM. LATITUDE, SMOOTHED SOLAR FLUX AND CRITICAL FREQUENCY
C (SEE MAIN PROGRAM).
C [REF.:K.RAWER,S.RAMAKRISHNAN,1978]
c NeQuick formula (Coisson et al, 2006), modified for computational
c speed and elegance (RS).
c----------------------------------------------------------------
      COMMON  /BLOCK1/HMF2,XNMF2,HMF1,F1REG
      common  /BLO11/B2TOP,TC3,hcor1
      logical 	f1reg

      DH=(H-HMF2)/B2TOP
      Z=DH/(1.0+DH/(8.0+0.01*DH))
      IF (Z.GT.40.0) THEN
         XE1=0.0
	 RETURN
      ENDIF
      Z=EXP(Z)
      IF (Z.GT.1E7) THEN
         Z=4.0/Z
      ELSE
         Z=4.0*Z/(1.0+Z)**2
      ENDIF
      XE1=XNMF2*Z
      RETURN
      END
C
C
      REAL FUNCTION XE2(H)
C ELECTRON DENSITY FOR THE BOTTOMSIDE F-REGION (HMF1...HMF2).
      COMMON    /BLOCK1/HMF2,XNMF2,HMF1,F1REG
      common    /BLOCK2/B0,B1,C1
      common	/ARGEXP/ARGMAX
      logical	f1reg

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
C ELECTRON DENSITY FOR THE F1-LAYER (HZ.....HMF1)
C USING THE NEW DEFINED F1-LAYER FUNCTION (Reinisch and Huang, Advances
C in Space Research, Volume 25, Number 1, 81-88, 2000)
      COMMON	/BLOCK1/	HMF2,XNMF2,HMF1,F1REG
      common	/BLOCK2/	B0,B1,D1F1
      logical	f1reg
C
      h1bar=h
      if (f1reg) H1BAR=HMF1*(1.0-((HMF1-H)/HMF1)**(1.0+D1F1))
      XE3=XE2(H1BAR)
      RETURN
      END
C
C
      REAL FUNCTION XE4(H)
C ELECTRON DENSITY FOR THE INTERMEDIATE REGION (HEF...HZ)
C USING THE NEW DEFINED FUNCTION
      COMMON	/BLOCK3/	HZ,T,HST
      common	/BLOCK4/	HME,XNME,HEF
C
      if(hst.lt.0.0) then
	 XE4=xnme+t*(h-hef)
	 return
      endif
      IF(HST.EQ.HEF) THEN
         H1BAR=H
      ELSE
         H1BAR=HZ+0.5*T-SIGN(1.0,T)*SQRT(T*(0.25*T+HZ-H))
      ENDIF
      XE4=XE3(H1BAR)
      RETURN
      END
C
C
      REAL FUNCTION XE5(H)
C ELECTRON DENSITY FOR THE E AND VALLEY REGION (HME..HEF).
      LOGICAL NIGHT
      COMMON    /BLOCK4/        HME,XNME,HEF
      common	/BLOCK5/        NIGHT,E(4)
      T3=H-HME
      T1=T3*T3*(E(1)+T3*(E(2)+T3*(E(3)+T3*E(4))))
      IF (NIGHT) then
         XE5=XNME*EXP(T1)
      else
         XE5=XNME*(1+T1)
      endif
      RETURN
      END
C
C
      REAL FUNCTION XE6(H)
C ELECTRON DENSITY FOR THE D REGION (HA...HME).
      COMMON    /BLOCK4/        HME,XNME,HEF
      common    /BLOCK6/        HMD,XNMD,HDX
      common    /BLOCK7/        D1,XKK,FP30,FP3U,FP1,FP2
      IF (H.LE.HDX) THEN
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
      COMMON    /BLOCK1/HMF2,XNMF2,XHMF1,F1REG
      common    /BLOCK3/HZ,T,HST
      common    /BLOCK4/HME,XNME,HEF
      logical 	f1reg
      IF (F1REG) THEN
         HMF1=XHMF1
      ELSE
         HMF1=HMF2
      ENDIF
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
      ENDIf
      RETURN
      END

C*************************************************************
C************* PEAK VALUES ELECTRON DENSITY ******************
C*************************************************************
C
C
      real FUNCTION FOUT(XMODIP,XLATI,XLONGI,UT,FF0)
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
      real FUNCTION XMOUT(XMODIP,XLATI,XLONGI,UT,XM0)
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
      fof1ed=0.0
      if (chi.gt.90.0) return

      DLA =  YLATI
      F0 = 4.35 + DLA * ( 0.0058 - 1.2E-4 * DLA )
      F100 = 5.348 + DLA * ( 0.011 - 2.3E-4 * DLA )
      FS = F0 + ( F100 - F0 ) * R / 100.0
      XMUE = 0.093 + DLA * ( 0.0046 - 5.4E-5 * DLA ) + 3.0E-4 * R
      FOF1 = FS * COS( CHI * UMR ) ** XMUE
      CHI0 = 49.84733 + 0.349504 * DLA
      CHI100 = 38.96113 + 0.509932 * DLA
      CHIM = ( CHI0 + ( CHI100 - CHI0 ) * R / 100. )
      IF(CHI.GT.CHIM) FOF1=-FOF1
      FOF1ED = FOF1
      RETURN
      END
C
C
      real FUNCTION f1_c1(xmodip,hour,suxnon,saxnon)
c F1 layer shape parameter C1 after Reinisch and Huang, Advances in
c Space Research, Volume 25, Number 1, 81-88, 2000.

      common	/const/umr
      pi = umr * 180.

      ABSMDP=ABS(XMODIP)
      DELA=4.32
      IF(ABSMDP.GE.18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)

      C1OLD = 0.09 + 0.11/DELA
      if(suxnon.eq.saxnon) then
          c1 = 2.5 * c1old
      else
          c1 = 2.5*c1old*cos((HOUR-12.)/(suxnon-saxnon)*pi)
      endif
      if(c1.lt.0.0) c1=0.0
      f1_c1=c1
      return
      end
c
c
      SUBROUTINE f1_prob (sza,glat,rz12,f1prob,f1probl)
c Occurrence probability of F1 layer after Scotto et al., Advances in
c Space Research, Volume 20, Number 9, 1773-1775, 1997.
c Input: solar zenith angle (sza) in degrees, geomagnetic latitude
c (glat) in degrees, 12-month running mean of sunspot number (Rz12).
c Output: F1 occurrence probability without L-condition cases (f1prob)
c and with L-condition cases (f1probl)
      common /const/umr

      xarg = 0.5 + 0.5 * cos(sza*umr)
      a = 2.98 + 0.0854 * rz12
      b = 0.0107 - 0.0022 * rz12
      c = -0.000256 + 0.0000147 * rz12
      gamma = a + ( b + c * glat) * glat
      f1pr = xarg ** gamma
      if(f1pr.lt.1.e-3) f1pr=0.0
      f1prob=f1pr
      f1prl = xarg ** 2.36
      if(f1prl.lt.1.e-3) f1prl=0.0
      f1probl=f1prl
      return
      end
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
c
      if(xhi.ge.90) goto 100
      Y = 6.05E8 + 0.088E8 * R
      yy = cos ( xhi * umr )
      ymd = y * exp( -0.1 / ( yy**2.7 ) )
      if (ymd.lt.yw) ymd = yw
      xmded=ymd
      RETURN

100   XMDED=YW
      RETURN
      END
C
C
      REAL FUNCTION GAMMA1(SMODIP,SLAT,SLONG,HOUR,
     &                        IHARM,NQ,K1,M,MM,M3,SFE)
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
      DO 250 I=2,IHARM
      C(I)=C(1)*C(I-1)-S(1)*S(I-1)
      S(I)=C(1)*S(I-1)+S(1)*C(I-1)
250   CONTINUE
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
      DO 350 J=1,INDEX
      SUM=SUM+COEF(1+J)*SS
      XSINX(J+1)=SS
      SS=SS*S3
350   CONTINUE
      XSINX(NQ(1)+2)=SS
      NP=NQ(1)+1
      SS=COS(SLAT*UMR)
      S3=SS
      DO 400 J=2,K1
      S0=SLONG*(J-1.)*UMR
      S1=COS(S0)
      S2=SIN(S0)
      INDEX=NQ(J)+1
      DO 450 L=1,INDEX
      NP=NP+1
      SUM=SUM+COEF(NP)*XSINX(L)*SS*S1
      NP=NP+1
      SUM=SUM+COEF(NP)*XSINX(L)*SS*S2
450   CONTINUE
      SS=SS*S3
400   CONTINUE
      GAMMA1=SUM
      RETURN
      END
C
C
C************************************************************
C***************** PROFILE PARAMETERS ***********************
C************************************************************
C
C
      REAL FUNCTION B0_98 ( HOUR, SAX, SUX, SEASON, R, ZLO, ZMODIP)
C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON is northern season with
C ISEASON=1 northern spring), low and high solar activity Rz12=10,
C 100 (IR=1,2), and modified dip latitudes of 0, 18 and 45
C degress (ILATI=1,2,3). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C
C input:
C       hour    LT in decimal hours
C       SAX     time of sunrise in decimal hours
C       SUX     time of sunset in decimal hours
C       SEASON  season in northern hemisphere (1=spring)
C               (real number between 0 and 4)
C       R       12-month running mean of sunspot number
C       ZLO     longitude
C       ZMODIP  modified dip latitude
C
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
C Updates (B0_new -> B0_98):
C
C 01/98 corrected to include a smooth transition at the modip equator
C       and no discontinuity at the equatorial change in season.
C 09/98 new B0 values incl values at the magnetic equator
C 10/98 longitude as input to determine if magnetic equator in northern
C         or southern hemisphere
C
      REAL	NITVAL
      DIMENSION B0F(2,4,2,3),bfd(2,3),zx(5),g(6),dd(5)
      DATA      B0F/201,68,210,61,192,68,199,67,240,80,245,83,
     &              233,71,230,65,108,65,142,81,110,68,77,75,
     &              124,98,164,100,120,94,96,112,78,81,94,84,
     &              81,81,65,70,102,87,127,91,109,88,81,78/
      DIMENSION BFR(2,4)
      data    zx/45.,72.,90.,108.,135./,dd/5*3.0/

!#ifdef B0_ORIG
!      nseasn=int(season)
!#else
      seas=season+3.5
      nseasn=int(seas)
      seas=seas-nseasn
!#endif
      zz = zmodip + 90.
      zz0 = 0.

C Interpolation in Rz12: linear from 10 to 100
      DO ISL=1,3
	  do isd=1,2
	     do iseas=1,4
                bfr(isd,iseas) = b0f(isd,iseas,1,isl) +
     &          (b0f(isd,iseas,2,isl)-b0f(isd,iseas,1,isl))/90.*(R-10.)
	     enddo
	  enddo
C Interpolation day/night with transitions at SAX (sunrise)
C and SUX (sunset) for northern/southern hemisphere iss=1/2
	  do iss=1,2
	     iseas = mod(nseasn + iss*2 + 1,4) + 1
!#ifdef B0_ORIG
!	     DAYVAL = bfr(1,iseas)
!	     NITVAL = bfr(2,iseas)
!#else
	     jseas = mod(nseasn + iss*2 + 2,4) + 1
	     DAYVAL = bfr(1,iseas) * (1.0 - seas) + bfr(1,jseas) * seas
	     NITVAL = bfr(2,iseas) * (1.0 - seas) + bfr(2,jseas) * seas
!#endif
             BFD(iss,ISL) = HPOL(HOUR,DAYVAL,NITVAL,SAX,SUX,1.,1.)
	  enddo
      enddo

C Interpolation with epstein-transitions in modified dip latitude.
C Transitions at +/-18 and +/-45 degrees; constant above +/-45.
C
C g(1:5) are the latitudinal slopes of B0;
C       g(1) is for the region from -90 to -45 degrees
C       g(2) is for the region from -45 to -18 degrees
C       g(3) is for the region from -18 to   0 degrees
C       g(4) is for the region from   0 to  18 degrees
C       g(5) is for the region from  18 to  45 degrees
C       g(6) is for the region from  45 to  90 degrees
C
C B0 =  bfd(2,3) at modip = -45,
C       bfd(2,2) at modip = -18,
C       bfd(2,1) or bfd(1,1) at modip = 0,
C       bfd(1,2) at modip = 20,
C       bfd(1,3) at modip = 45.
C
C If the Longitude is between 200 and 320 degrees than the modip
C equator is in the southern hemisphere and bfd(2,1) is used at the
C equator, otherwise bfd(1,1) is used.
c
      zx1=bfd(2,3)
      zx2=bfd(2,2)
      zx3=bfd(1,1)
      if(zlo.gt.200.0.and.zlo.lt.320.0) zx3=bfd(2,1)
      zx4=bfd(1,2)
      zx5=bfd(1,3)
      g(1) = 0.
      g(2) = ( zx2 - zx1 ) / 27.
      g(3) = ( zx3 - zx2 ) / 18.
      g(4) = ( zx4 - zx3 ) / 18.
      g(5) = ( zx5 - zx4 ) / 27.
      g(6) = 0.
      B0_98 = zx1
      DO I=1,5
        aa = eptr(zz ,dd(i),zx(i))
        bb = eptr(zz0,dd(i),zx(i))
        DSUM = (G(I+1) - G(I)) * (AA-BB) * dd(i)
        B0_98 = B0_98 + DSUM
      enddo

      RETURN
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
      LOGICAL L1/.FALSE./,LINKS,K,SCHALT
      SCHALT=.FALSE.
      EP=EPS
      X1=X11
      X2=X22
      F1=FX11-FW
      F2=FX22-FW
      K=.FALSE.
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
      IF((LINKS.AND.(.NOT.L1)).OR.(.NOT.LINKS.AND.L1)) NG=2*NG
      GOTO 200
800   RETURN
      END
C
C
C******************************************************************
C********** ZENITH ANGLE, DAY OF YEAR, TIME ***********************
C******************************************************************
C
C
      SUBROUTINE SOCO (ld,t,flat,Elon,height,
     &          DECLIN, ZENITH, SUNRSE, SUNSET)
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
c		height	height in km
c
c out:  declin      declination of the sun in degrees
c       zenith      zenith angle of the sun in degrees
c       sunrse      local time of sunrise in hours
c       sunset      local time of sunset in hours
c-------------------------------------------------------------------
c
      common /const/   umr
c amplitudes of Fourier coefficients  --  1955 epoch.................
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
      phi = 15. * umr * ( t - 12.) + et
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
      h=height*1000.
      chih = 90.83 + 0.0347 * sqrt(h)
c this includes corrections for horizontal refraction and........
c semi-diameter of the solar disk................................
      ch = cos(chih * umr)
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
      et = et / umr / 15.
      phi = phi / umr / 15.
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
      IF(ABS(SU).GT.25.) THEN
         IF(SU.GT.0.0) THEN
            HPOL=TW
         ELSE
            HPOL=XNW
         ENDIF
         RETURN
      ENDIF
      HPOL=XNW+(TW-XNW)*EPST(HOUR,DSA,SA)+
     &  (XNW-TW)*EPST(HOUR,DSU,SU)
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
c  years also evenly divisible by 100 are not leap years, except years
c  also evenly divisible by 400 are leap years. The year 2000 therefore
C  is a leap year. The 100 and 400 year exception rule
c     if((iyear/4*4.eq.iyear).and.(iyear/100*100.ne.iyear)) mm(2)=29
c  will become important again in the year 2100 which is not a leap
C  year.
c
      mm(2)=28
      if(iyear/4*4.eq.iyear) mm(2)=29

      IF (IN.LE.0) THEN
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
C *********************************************************************
C ************************ EPSTEIN FUNCTIONS **************************
C *********************************************************************
C REF:  H. G. BOOKER, J. ATMOS. TERR. PHYS. 39, 619-623, 1977
C       K. RAWER, ADV. SPACE RES. 4, #1, 11-15, 1984
C *********************************************************************
C
C
      REAL FUNCTION EPTR ( X, SC, HX )
C --------------------------------------------------------- TRANSITION
      COMMON/ARGEXP/ARGMAX
      D1 = ( X - HX ) / SC
      IF (ABS(D1).LT.ARGMAX) GOTO 1
      IF (D1.GT.0.0) THEN
        EPTR = D1
      ELSE
        EPTR = 0.0
      ENDIF
      RETURN
1     EPTR = ALOG ( 1. + EXP( D1 ))
      RETURN
      END
C
C
      REAL FUNCTION EPST ( X, SC, HX )
C -------------------------------------------------------------- STEP
      COMMON/ARGEXP/ARGMAX
      D1 = ( X - HX ) / SC
      IF (ABS(D1).LT.ARGMAX) GOTO 1
      IF (D1.GT.0.0) THEN
        EPST = 1.
      ELSE
        EPST = 0.
      ENDIF
      RETURN
1     EPST = 1. / ( 1. + EXP( -D1 ))
      RETURN
      END
C
C
      SUBROUTINE tcon(yr,mm,day,idn,rz,ig,rsn,nmonth)
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
c the indices are obtained from the indices file ig_rz.dat that is
c read in subroutine initialize and stored in COMMON/indices/
c----------------------------------------------------------------

      integer      yr, mm, day, iflag, iyst, iyend,iymst
      integer      imst,iymend
      integer	   unit,freeunit
      character*80 dirnam
      real         ionoindx(1202),indrz(1202)
      real         ig(3),rz(3)

      common /iounit/konsol,kerror

      save         ionoindx,indrz,iflag,iyst,iymst,
     &             iymend,imst
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
c month following the end month. These values are also included in
c IG_RZ.
c
c A negative Rz index means that the given index is the 13-months-
C running mean of the solar radio flux (F10.7). The close correlation
C between (Rz)12 and (F10.7)12 is used to derive the (Rz)12 indices.
c
c An IG index of -111 indicates that no IG values are available for the
c time period. In this case a correlation function between (IG)12 and
C (Rz)12 is used to obtain (IG)12.
c
c The computation of the 13-month-running mean for month M requires the
c indices for the six months preceeding M and the six months following
C M (month: M-6, ..., M+6). To calculate the current running mean one
C therefore requires predictions of the indix for the next six months.
C Starting from six months before the UPDATE DATE (listed at the top of
c the file) and onward the indices are therefore based on indices
c predictions.
c
      if(iflag.eq.0) then

          unit=freeunit()
          dirnam='/user/altim'
          call checkenv('ALTIM',dirnam,l)
	  if (konsol.ge.0) write (konsol,*) "Reading: ig_rz.dat"
          open(unit,file=dirnam(:l)//'/data/indices/ig_rz.dat',
     |		status='old')

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
          do 1 jj=1,inum_vals
                rrr=indrz(jj)
                if(rrr.lt.0.0) then
                        covr=abs(rrr)
                        rrr=33.52*sqrt(covr+85.12)-408.99
                        if(rrr.lt.0.0) rrr=0.0
                        indrz(jj)=rrr
                        endif
                if(ionoindx(jj).gt.-90.) goto 1
                  zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
                  if(zi.gt.274.0) zi=274.0
                  ionoindx(jj)=zi
1               continue
          close(unit)
          iflag = 1
        endif

        iytmp=yr*100+mm
        if (iytmp .lt. iymst .or. iytmp .gt. iymend) then
               if(konsol.ge.0) write(konsol,8000) iytmp,iymst,
     &                                            iymend
 8000          format(1x,I10,'** OUT OF RANGE **'/,5x,
     &  'The file IG_RZ.DAT which contains the indices Rz12',
     &  ' and IG12'/5x,'currently only covers the time period',
     &  ' (yymm) : ',I6,'-',I6)
               nmonth=-1
               return
               endif

        num=2-imst+(yr-iyst)*12+mm

        rz(1)=indrz(num)
        ig(1)=ionoindx(num)
        midm=15
        if(mm.eq.2) midm=14
        call MODA(0,yr,mm,midm,idd1,nrdaym)
        if(day.lt.midm) goto 1926
                imm2=mm+1
                if(imm2.gt.12) then
                        imm2=1
                        iyy2=yr+1
                        idd2=380
                        if(yr/4*4.eq.yr) idd2=381

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
                goto 1927
1926            imm2=mm-1
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

1927    nmonth=imm2
            return
            end
c
c
        SUBROUTINE APF(IYYYY,IMN,ID,HOUR,IAP)
c--------------------------------------------------------------------
c Finds 3-hourly Ap indices for IRI-storm for given year IYYYY (yyyy),
c month (IMN), day (ID), and UT (HOUR, decimal hours). The indices are
c stored in IAP(13) providing the 13 3-hourly indices prior to HOUR.
c The 3-hour UT intervals during the day are: (0-3),)3-6),)6-9),9-12,
c 12-15,15-18,18-21,)21-24(.
c If date is outside the range of the Ap indices file than iap(1)=-5
c--------------------------------------------------------------------

      parameter	(map=200000)
      DIMENSION iiap(map),iap(13),lm(12)
      integer unit,freeunit,nap/0/,i
      character*80 dirnam
      common /iounit/konsol,kerror
      save nap
      DATA LM/31,28,31,30,31,30,31,31,30,31,30,31/

      if (nap.eq.0) then
         unit=freeunit()
         dirnam='/user/altim'
         call checkenv('ALTIM',dirnam,l)
	 if (konsol.ge.0) write (konsol,*) "Reading: ap_flux.dat"
         Open(unit,FILE=dirnam(:l)//'/data/indices/ap_flux.dat',
     *  STATUS='OLD')
20	 read (unit,'(8x,8i3)',end=90) (iiap(nap+i),i=1,8)
	 nap=nap+8
	 goto 20
90	 close(unit)
      endif

      do i=1,8
         iap(i)=-1
      enddo

      if(iyyyy.lt.1958) goto 21   ! file starts at Jan 1, 1958
      is=0
      do i=1958,iyyyy-1
         nyd=365
         if(i/4*4.eq.i) nyd=366
         IS=IS+nyd
      enddo

      LM(2)=28
      if(iyyyy/4*4.eq.iyyyy) LM(2)=29
      do i=1,IMN-1
         IS=IS+LM(i)
      ENDDO
      IS=IS+ID-1

      ihour=int(hour/3.)+1
      if(ihour.gt.8) ihour=8
      is=is*8+ihour

      if(is.lt.13 .or. is.gt.nap) goto 21   ! at least 13 indices available

      do i=1,13
         iap(i)=iiap(is-13+i)
      enddo
      return

21    if (konsol.ge.0) write(konsol,100)
100   format(1X,'Date is outside range of Ap indices file.',
     &   ' STORM model is turned off.')
      IAP(1)=-5
      RETURN
      END
C
C
C----------------------STORM MODEL --------------------------------
C
      SUBROUTINE CONVER(rga,rgo,rgma)

C     This subroutine converts a geographic latitude and longitude
C     location to a corrected geomagnetic latitude.
C
C     INPUT:
C       geographic latitude   -90. to +90.
C       geographic longitude  0. to 360. positive east from Greenwich.
C
C     OUTPUT:
C       corrected geomagnetic latitude	-90. to +90.

      DIMENSION CORMAG(20,91)
      DATA ((CORMAG(i,j),i=1,20),j=1,31)/
     +163.68,163.68,163.68,163.68,163.68,163.68,
     +163.68,163.68,163.68,163.68,163.68,163.68,163.68,163.68,
     +163.68,163.68,163.68,163.68,163.68,163.68,162.60,163.12,
     +163.64,164.18,164.54,164.90,165.16,165.66,166.00,165.86,
     +165.20,164.38,163.66,162.94,162.42,162.00,161.70,161.70,
     +161.80,162.14,161.20,162.18,163.26,164.44,165.62,166.60,
     +167.42,167.80,167.38,166.82,166.00,164.66,163.26,162.16,
     +161.18,160.40,159.94,159.80,159.98,160.44,159.80,161.14,
     +162.70,164.50,166.26,167.90,169.18,169.72,169.36,168.24,
     +166.70,164.80,162.90,161.18,159.74,158.60,157.94,157.80,
     +157.98,158.72,158.40,160.10,162.02,164.28,166.64,169.00,
     +170.80,171.72,171.06,169.46,167.10,164.64,162.18,160.02,
     +158.20,156.80,156.04,155.80,156.16,157.02,157.00,158.96,
     +161.24,163.86,166.72,169.80,172.42,173.72,172.82,170.34,
     +167.30,164.22,161.34,158.74,156.60,155.00,154.08,153.90,
     +154.36,155.36,155.50,157.72,160.36,163.32,166.60,170.20,
     +173.70,175.64,174.18,170.80,167.10,163.56,160.24,157.36,
     +154.96,153.10,152.08,151.92,152.46,153.76,154.10,156.52,
     +159.36,162.52,166.24,170.30,174.62,177.48,175.04,170.82,
     +166.60,162.70,159.02,155.88,153.22,151.20,150.08,149.92,
     +150.64,152.20,152.80,155.32,158.28,161.70,165.58,170.00,
     +174.84,178.46,175.18,170.38,165.80,161.64,157.80,154.38,
     +151.52,149.30,148.18,148.02,148.92,150.60,151.40,154.08,
     +157.18,160.68,164.78,169.40,174.34,177.44,174.28,169.44,
     +164.70,160.34,156.30,152.78,149.72,147.40,146.18,146.04,
     +147.12,149.04,150.10,152.88,156.00,159.58,163.78,168.50,
     +173.28,175.60,172.86,168.14,163.40,158.98,154.88,151.10,
     +147.98,145.50,144.18,144.14,145.40,147.48,148.80,151.68,
     +154.88,158.48,162.68,167.40,171.76,173.60,171.12,166.68,
     +162.00,157.48,153.28,149.50,146.18,143.50,142.18,142.24,
     +143.68,145.98,147.50,150.54,153.68,157.28,161.42,166.10,
     +170.10,171.48,169.22,164.98,160.40,155.88,151.68,147.80,
     +144.34,141.60,140.18,140.26,141.98,144.62,146.30,149.34,
     +152.48,155.98,160.08,164.60,168.34,169.38,167.20,163.18,
     +158.60,154.18,149.98,146.02,142.54,139.70,138.18,138.46,
     +140.26,143.16,145.10,148.14,151.18,154.60,158.68,163.10,
     +166.48,167.28,165.18,161.32,156.90,152.48,148.28,144.32,
     +140.74,137.80,136.22,136.48,138.64,141.76,143.90,146.98,
     +149.98,153.30,157.24,161.40,164.52,165.16,162.86,159.42,
     +155.00,150.68,146.48,142.52,138.94,135.90,134.22,134.68,
     +137.02,140.40,142.70,145.84,148.76,151.92,155.74,159.70,
     +162.52,162.96,160.98,157.42,153.10,148.84,144.68,140.82,
     +137.20,134.00,132.32,132.80,135.42,139.10,141.60,144.74,
     +147.46,150.52,154.20,158.00,160.46,160.76,158.86,155.36,
     +151.20,146.94,142.88,139.02,135.40,132.10,130.32,131.00,
     +133.80,137.74,140.50,143.58,146.24,149.12,152.60,156.20,
     +158.40,158.66,156.76,153.36,149.30,145.04,141.08,137.30,
     +133.60,130.30,128.42,129.12,132.28,136.44,139.30,142.48,
     +144.94,147.64,150.48,154.30,156.34,156.36,154.56,151.26,
     +147.30,143.14,139.20,135.50,131.90,128.40,126.52,127.32,
     +130.76,135.18,138.20,141.28,143.72,146.24,149.26,152.40,
     +154.24,154.16,152.36,149.16,145.30,141.24,137.30,133.70,
     +130.10,126.60,124.62,125.54,129.16,133.92,137.10,140.18,
     +142.42,144.66,147.62,150.50,152.18,151.96,150.16,147.10,
     +143.30,139.24,135.50,131.90,128.36,124.80,122.72,123.74,
     +127.64,132.62,135.90,139.02,141.12,143.18,145.92,148.60,
     +149.98,149.76,148.04,145.00,141.20,137.30,133.60,130.10,
     +126.60,123.00,120.86,121.96,126.12,131.36,134.80,137.88,
     +139.80,141.68,144.08,146.60,147.88,147.56,145.84,142.90,
     +139.20,135.30,131.70,128.28,124.86,121.30,118.96,120.18,
     +124.70,130.16,133.60,136.72,138.48,140.10,142.38,144.60,
     +145.72,145.34,143.64,140.80,137.10,133.30,129.72,126.48,
     +123.10,119.50,117.16,118.48,123.18,128.86,132.40,135.42,
     +137.08,138.50,140.54,142.60,143.52,143.06,141.44,138.70,
     +135.10,131.30,127.82,124.58,121.40,117.70,115.26,116.70,
     +121.66,127.60,131.20,134.22,135.66,136.82,138.70,140.60,
     +141.36,140.86,139.24,136.50,133.00,129.30,125.92,122.78,
     +119.60,116.00,113.40,114.92,120.16,126.30,130.00,132.92,
     +134.24,135.14,136.80,138.60,139.16,138.64,137.12,134.40,
     +130.90,127.20,123.92,120.96,117.90,114.20,111.56,113.12,
     +118.64,124.90,128.70,131.56,132.74,133.44,134.90,136.50,
     +137.00,136.36,134.82,132.30,128.70,125.16,121.94,119.06,
     +116.10,112.50,109.70,111.42,117.14,123.60,127.30,130.16,
     +131.22,131.66,133.00,134.50,134.80,134.14,132.62,130.14,
     +126.60,123.06,119.94,117.16,114.30,110.70,107.80,109.64,
     +115.62,122.24,125.90,128.76,129.62,129.96,131.06,132.40,
     +132.60,131.86,130.42,128.00,124.50,120.96,117.96,115.26,
     +112.54,108.90,105.94,107.86,114.02,120.84/

      DATA ((CORMAG(i,j),i=1,20),j=32,61)/
     +124.05,126.79,
     +127.55,127.83,128.90,130.21,130.41,129.71,128.33,125.96,
     +122.49,118.96,115.97,113.26,110.52,106.89,104.01,106.00,
     +112.21,119.06,122.19,124.82,125.48,125.69,126.73,128.03,
     +128.22,127.55,126.23,123.92,120.47,116.97,113.97,111.26,
     +108.50,104.89,102.08,104.14,110.41,117.29,120.34,122.85,
     +123.40,123.56,124.57,125.84,126.03,125.40,124.14,121.88,
     +118.46,114.97,111.98,109.26,106.48,102.88,100.15,102.28,
     +108.60,115.51,118.49,120.88,121.33,121.42,122.40,123.65,
     +123.84,123.24,122.04,119.83,116.45,112.97,109.98,107.26,
     +104.46,100.87,098.22,100.42,106.79,113.74,116.63,118.91,
     +119.26,119.29,120.24,121.47,121.65,121.09,119.95,117.79,
     +114.43,110.98,107.99,105.26,102.44,098.87,096.29,098.56,
     +104.98,111.96,114.78,116.94,117.19,117.15,118.07,119.28,
     +119.46,118.93,117.86,115.75,112.42,108.98,106.00,103.26,
     +100.42,096.86,094.36,096.70,103.18,110.19,112.93,114.97,
     +115.12,115.02,115.91,117.09,117.27,116.78,115.76,113.71,
     +110.41,106.98,104.00,101.26,098.40,094.85,092.43,094.84,
     +101.37,108.41,111.07,113.00,113.04,112.88,113.74,114.91,
     +115.08,114.62,113.67,111.67,108.39,104.99,102.01,099.26,
     +096.38,092.85,090.51,092.97,099.56,106.64,109.22,111.03,
     +110.97,110.75,111.58,112.72,112.89,112.47,111.57,109.63,
     +106.38,102.99,100.01,097.26,094.36,090.84,088.58,091.11,
     +097.75,104.86,107.37,109.06,108.90,108.61,109.41,110.53,
     +110.70,110.31,109.48,107.59,104.37,100.99,098.02,095.26,
     +092.34,088.83,086.65,089.25,095.95,103.09,105.51,107.09,
     +106.83,106.48,107.25,108.35,108.51,108.16,107.39,105.55,
     +102.35,099.00,096.03,093.26,090.32,086.83,084.72,087.39,
     +094.14,101.31,103.66,105.12,104.76,104.34,105.08,106.16,
     +106.32,106.00,105.29,103.50,100.34,097.00,094.03,091.26,
     +088.30,084.82,082.79,085.53,092.33,099.54,101.81,103.15,
     +102.68,102.21,102.92,103.97,104.13,103.85,103.20,101.46,
     +098.33,095.00,092.04,089.26,086.28,082.81,080.86,083.67,
     +090.52,097.76,099.95,101.18,100.61,100.07,100.75,101.79,
     +101.94,101.69,101.10,099.42,096.31,093.01,090.04,087.26,
     +084.26,080.81,078.93,081.81,088.72,095.99,098.10,099.21,
     +098.54,097.94,098.59,099.60,099.75,099.54,099.01,097.38,
     +094.30,091.01,088.05,085.26,082.24,078.80,077.00,079.95,
     +086.91,094.21,096.25,097.24,096.47,095.81,096.43,097.41,
     +097.56,097.39,096.92,095.34,092.29,089.01,086.06,083.26,
     +080.22,076.79,075.07,078.09,085.10,092.43,094.39,095.27,
     +094.40,093.67,094.26,095.23,095.37,095.23,094.82,093.30,
     +090.27,087.02,084.06,081.26,078.20,074.79,073.14,076.23,
     +083.30,090.66,092.54,093.30,092.32,091.54,092.10,093.04,
     +093.18,093.08,092.73,091.26,088.26,085.02,082.07,079.26,
     +076.18,072.78,071.21,074.37,081.49,088.88,090.69,091.33,
     +090.25,089.40,089.93,090.85,090.99,090.92,090.63,089.21,
     +086.25,083.02,080.07,077.26,074.16,070.77,069.28,072.51,
     +079.68,087.11,088.83,089.36,088.18,087.27,087.77,088.67,
     +088.80,088.77,088.54,087.17,084.23,081.03,078.08,075.26,
     +072.14,068.77,067.35,070.65,077.87,085.33,086.98,087.39,
     +086.11,085.13,085.60,086.48,086.61,086.61,086.45,085.13,
     +082.22,079.03,076.09,073.26,070.12,066.76,065.42,068.79,
     +076.07,083.56,085.13,085.42,084.04,083.00,083.44,084.29,
     +084.42,084.46,084.35,083.09,080.21,077.03,074.09,071.26,
     +068.10,064.75,063.49,066.93,074.26,081.78,083.27,083.45,
     +081.96,080.86,081.27,082.11,082.23,082.30,082.26,081.05,
     +078.19,075.04,072.10,069.26,066.08,062.75,061.57,065.06,
     +072.45,080.01,081.42,081.48,079.89,078.73,079.11,079.92,
     +080.04,080.15,080.16,079.01,076.18,073.04,070.10,067.26,
     +064.06,060.74,059.64,063.20,070.64,078.23,079.57,079.51,
     +077.82,076.59,076.94,077.73,077.85,077.99,078.07,076.97,
     +074.17,071.04,068.11,065.26,062.04,058.73,057.71,061.34,
     +068.84,076.46,077.71,077.54,075.75,074.46,074.78,075.55,
     +075.66,075.84,075.98,074.93,072.15,069.05,066.12,063.26,
     +060.02,056.73,055.78,059.48,067.03,074.68,075.86,075.57,
     +073.68,072.32,072.61,073.36,073.47,073.68,073.88,072.88,
     +070.14,067.05,064.12,061.26,058.00,054.72,053.85,057.62,
     +065.22,072.91,074.01,073.60,071.60,070.19,070.45,071.17,
     +071.28,071.53,071.79,070.84,068.13,065.05,062.13,059.26,
     +055.98,052.71,051.92,055.76,063.41,071.13,072.15,071.63,
     +069.53,068.05,068.28,068.99,069.09,069.37,069.69,068.80,
     +066.11,063.06,060.13,057.26,053.96,050.71,049.99,053.90,
     +061.61,069.36,070.30,069.66,067.46,065.92,066.12,066.80,
     +066.90,067.22,067.60,066.76,064.10,061.06,058.14,055.26,
     +051.94,048.70,048.06,052.04,059.80,067.58/

      DATA ((CORMAG(i,j),i=1,20),j=62,91)/
     +067.70,067.06,
     +065.08,063.72,063.98,064.60,064.80,065.12,065.60,064.86,
     +062.40,059.26,056.24,053.18,049.84,046.60,046.12,050.12,
     +057.52,064.80,064.90,064.42,062.70,061.62,061.78,062.40,
     +062.60,063.04,063.58,063.00,060.60,057.46,054.42,051.18,
     +047.70,044.60,044.22,048.02,055.06,061.92,062.10,061.72,
     +060.32,059.50,059.68,060.20,060.46,060.94,061.58,061.00,
     +058.70,055.66,052.52,049.18,045.60,042.50,042.22,046.00,
     +052.60,058.98,059.20,059.18,058.12,057.32,057.48,058.00,
     +058.30,058.84,059.48,059.04,056.90,053.86,050.62,047.10,
     +043.50,040.50,040.28,043.98,050.22,056.18,056.40,056.64,
     +055.84,055.20,055.38,055.80,056.16,056.84,057.48,057.04,
     +055.10,052.06,048.70,045.10,041.40,038.40,038.28,041.88,
     +047.94,053.44,053.70,054.14,053.56,053.10,053.24,053.70,
     +054.06,054.74,055.38,055.14,053.20,050.26,046.80,043.10,
     +039.34,036.40,036.38,039.96,045.56,050.84,051.10,051.70,
     +051.36,051.00,051.14,051.50,051.96,052.64,053.38,053.08,
     +051.30,048.36,044.90,041.02,037.24,034.40,034.38,037.86,
     +043.28,048.20,048.50,049.26,049.18,048.90,049.04,049.40,
     +049.86,050.64,051.28,051.08,049.40,046.46,042.98,039.02,
     +035.14,032.40,032.48,035.72,041.00,045.70,046.00,046.96,
     +046.98,046.80,046.94,047.30,047.76,048.54,049.28,049.08,
     +047.40,044.56,041.08,037.02,033.14,030.40,030.58,033.84,
     +038.72,043.20,043.50,044.62,044.80,044.80,044.94,045.20,
     +045.76,046.54,047.18,046.98,045.50,042.66,039.08,035.02,
     +031.14,028.40,028.58,031.82,036.52,040.80,041.20,042.32,
     +042.54,042.70,042.84,043.20,043.66,044.44,045.08,044.98,
     +043.50,040.76,037.08,033.04,029.04,026.40,026.68,029.82,
     +034.34,038.40,038.80,040.12,040.60,040.70,040.84,041.10,
     +041.62,042.34,042.98,042.88,041.50,038.76,035.18,031.04,
     +027.14,024.50,024.78,027.70,032.14,036.06,036.50,037.88,
     +038.50,038.68,038.84,039.10,039.56,040.34,040.88,040.82,
     +039.40,036.76,033.18,029.12,025.14,022.50,022.88,025.90,
     +029.96,033.86,034.30,035.68,036.42,036.68,036.84,037.10,
     +037.56,038.24,038.88,038.72,037.40,034.76,031.18,027.12,
     +023.14,020.60,020.98,023.90,027.88,031.66,032.10,033.58,
     +034.32,034.68,034.84,035.10,035.56,036.24,036.78,036.62,
     +035.30,032.72,029.18,025.14,021.24,018.70,019.08,021.90,
     +025.88,029.42,029.90,031.48,032.32,032.68,032.84,033.10,
     +033.56,034.22,034.68,034.42,033.20,030.72,027.28,023.22,
     +019.34,016.80,017.24,020.00,023.78,027.32,027.70,029.38,
     +030.24,030.68,030.94,031.20,031.66,032.22,032.58,032.32,
     +031.10,028.62,025.28,021.32,017.48,015.00,015.38,018.18,
     +021.80,025.22,025.70,027.28,028.24,028.78,029.04,029.30,
     +029.66,030.22,030.50,030.22,029.00,026.62,023.30,019.42,
     +015.64,013.10,013.54,016.28,019.80,023.12,023.60,025.24,
     +026.24,026.78,027.14,027.40,027.76,028.22,028.40,028.12,
     +026.80,024.52,021.30,017.52,013.78,011.30,011.74,014.48,
     +017.90,021.12,021.60,023.24,024.34,024.88,025.24,025.50,
     +025.86,026.22,026.40,025.98,024.70,022.48,019.40,015.72,
     +012.04,009.50,009.94,012.58,016.02,019.12,019.60,021.24,
     +022.34,022.98,023.34,023.70,024.00,024.30,024.40,023.88,
     +022.60,020.48,017.52,014.00,010.34,007.80,008.18,010.88,
     +014.22,017.18,017.60,019.34,020.44,021.16,021.54,021.90,
     +022.16,022.40,022.32,021.78,020.60,018.48,015.62,012.20,
     +008.68,006.00,006.44,009.18,012.42,015.28,015.80,017.44,
     +018.54,019.26,019.74,020.10,020.30,020.50,020.32,019.72,
     +018.50,016.54,013.84,010.68,007.14,004.40,004.74,007.58,
     +010.74,013.48,014.00,015.54,016.74,017.46,017.94,018.30,
     +018.50,018.58,018.32,017.72,016.50,014.64,012.24,009.18,
     +005.84,002.90,003.30,006.16,009.14,011.84,012.30,013.78,
     +014.94,015.66,016.24,016.50,016.70,016.70,016.42,015.78,
     +014.60,012.90,010.66,007.86,004.88,001.60,001.72,004.96,
     +007.84,010.24,010.70,012.14,013.24,013.96,014.44,014.80,
     +014.90,014.88,014.52,013.92,012.80,011.30,009.28,006.94,
     +004.32,001.80,001.94,004.34,006.78,008.94,009.40,010.58,
     +011.64,012.36,012.74,013.10,013.20,013.08,012.72,012.12,
     +011.10,009.86,008.30,006.50,004.60,003.10,003.16,004.50,
     +006.20,007.90,008.40,009.42,010.14,010.76,011.14,011.40,
     +011.40,011.38,011.02,010.46,009.70,008.72,007.64,006.46,
     +005.42,004.60,004.70,005.34,006.24,007.36,007.90,008.46,
     +008.92,009.28,009.54,009.70,009.70,009.68,009.42,009.06,
     +008.60,008.08,007.56,007.02,006.56,006.30,006.30,006.52,
     +006.96,007.38,008.15,008.15,008.15,008.15,008.15,008.15,
     +008.15,008.15,008.15,008.15,008.15,008.15,008.15,008.15,
     +008.15,008.15,008.15,008.15,008.15,008.15/

C     Data Input
      rlan = rga
      rlo = rgo

C     From "normal" geographic latitude
C     to angle from South Pole.
      rla = rlan + 90

      IF (rlo .EQ. 360) THEN
      	rlo = 0
        END IF

C     PROXIMITY

C     coefficients of the latitudinal points
      LA1 = (INT(rla/2)+1)
      LA2 = LA1 + 1
      if(la2.gt.91) la2=91

C     coefficients of the longitudinal points
      LO1 = (INT(rlo/18)+1)
corr      LO2 = LO1 + 1
      LO2 = MOD(LO1,20) + 1

C     Four points of Geomagnetic Coordinates
      gm1 = CORMAG(LO1,LA1)
      gm2 = CORMAG(LO1,LA2)
      gm3 = CORMAG(LO2,LA1)
      gm4 = CORMAG(LO2,LA2)

C     latitudinal points
      X1 = ABS(rla - (INT(rla)))
      X2 = 2. - X1

C     longitudinal points
      Y1 = ABS(rlo - (INT(rlo)))
      Y2 = 18. - Y1

C     X AND Y VALUES
      x = X1 / (X1 + X2)
      y = Y1 / (Y1 + Y2)

C     INTERPOLATION
      gmla = gm1 * (1 - x) * (1 - y) + gm2 * (1 - y) * (x) + gm3 * (y)
     1 * (1 - x) + gm4 * (x) * (y)

C     OUTPUT OF THE PROGRAM
C     From corrected geomagnetic latitude from North Pole
C     to "normal"  geomagnetic latitude.
      rgma = 90. - gmla

      END
c
c
      SUBROUTINE STORM(ap,rga,rgo,coor,rgma,ut,doy,cf)
C----------------------------------------------------------------------
C      Fortran code to obtain the foF2 storm-time correction factor at
C      a given location and time, using the current and the 12 previous
C      ap values as input.
C
C      ap ---> (13 elements integer array). Array with the preceeding
C              13 value of the 3-hourly ap index. The 13th value
C              in the array will contain the ap at the UT of interest,
C              the 12th value will contain the 1st three hourly interval
C              preceeding the time of interest, and so on to the 1st
C              ap value at the earliest time.
C     coor --> (integer). If coor = 2, rga should contain the
C                         geomagnetic latitude.
C                         If coor = 1, rga should contain the
C                         geographic latitude.
C     rga ---> (real, -90 to 90) geographic or geomagnetic latitude.
C     rgo ---> (real, 0 to 360, positive east from Greenwich.)
C                           geographic longitude, only used if coor=1.
C     ut  ---> (integer, hours 00 to 23) Universal Time of interest.
C     doy ---> (integer, 1 to 366)Day of the year.
C     cf  ---> (real) The output; the storm-time correction factor used
C              to scale foF2, foF2 * cf.
C
C     This model and computer code was developed by E. Araujo-Pradere,
C     T. Fuller-Rowell and M. Condrescu, SEC, NOAA, Boulder, USA
C     Ref:
C     T. Fuller-Rowell, E. Araujo-Pradere, and M. Condrescu, An
C       Empirical Ionospheric Storm-Time Ionospheric Correction Model,
C       Adv. Space Res. 8, 8, 15-24, 2000.
C----------------------------------------------------------------------
C     DIMENSIONS AND COEFFICIENTS VALUES
      common /iounit/konsol,kerror

      DIMENSION c3(20)
      DATA c3/0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,-9.44E-12,
     +0.00E+00,3.04E-12,0.00E+00,9.32E-12,-1.07E-11,0.00E+00,0.00E+00,
     +0.00E+00,1.09E-11,0.00E+00,0.00E+00,0.00E+00,0.00E+00,-1.01E-11/

      DIMENSION c2(20)
      DATA c2/1.16E-08,0.00E+00,0.00E+00,-1.46E-08,0.00E+00,9.86E-08,
     +2.25E-08,-1.67E-08,-1.62E-08,-9.42E-08,1.17E-07,4.32E-08,3.97E-08,
     +3.13E-08,-8.04E-08,3.91E-08,2.58E-08,3.45E-08,4.76E-08,1.13E-07/

      DIMENSION c1(20)
      DATA c1/-9.17E-05,-1.37E-05,0.00E+00,7.14E-05,0.00E+00,-3.21E-04,
     +-1.66E-04,-4.10E-05,1.36E-04,2.29E-04,-3.89E-04,-3.08E-04,
     +-2.81E-04,-1.90E-04,4.76E-05,-2.80E-04,-2.07E-04,-2.91E-04,
     +-3.30E-04,-4.04E-04/

      DIMENSION c0(20)
      DATA c0/1.0136E+00,1.0478E+00,1.00E+00,1.0258E+00,1.00E+00,
     +1.077E+00,1.0543E+00,1.0103E+00,9.9927E-01,9.6876E-01,1.0971E+00,
     +1.0971E+00,1.0777E+00,1.1134E+00,1.0237E+00,1.0703E+00,1.0248E+00,
     +1.0945E+00,1.1622E+00,1.1393E+00/

      DIMENSION fap(36)
      DATA fap/0.,0.,0.037037037,0.074074074,0.111111111,0.148148148,
     |0.185185185,0.222222222,0.259259259,0.296296296,0.333333333,
     |0.37037037,0.407407407,0.444444444,0.481481481,0.518518519,
     |0.555555556,0.592592593,0.62962963,0.666666667,0.703703704,
     |0.740740741,0.777777778,0.814814815,0.851851852,0.888888889,
     |0.925925926,0.962962963,1.,0.66666667,0.33333334,0.,0.333333,
     |0.666666,1.,0.7/

      integer code(8,6)
      data code/3,4,5,4,3,2,1,2,3,2,1,2,3,4,5,4,8,7,6,7,8,9,10,9,
     *13,12,11,12,13,14,15,14,18,17,16,17,18,19,20,19,18,17,16,17,
     *18,19,20,19/

      INTEGER ape(39)
      INTEGER ap(13)
      INTEGER ut,doy,dayno,coor,s1,s2,l1,l2
      REAL rgma, rap, rga, rgo, rs, rl

C      CALLING THE PROGRAM TO CONVERT TO GEOMAGNETIC COORDINATES

      IF (coor .EQ. 1) THEN
         CALL CONVER (rga,rgo,rgma)
      ELSE IF (coor .EQ. 2) THEN
         rgma = rga
      ELSE
         WRITE (kerror,*)'Wrong Coordinates Selection -------- >>', coor
         return
      ENDIF

      if(doy.gt.366.or.doy.lt.1)then
          WRITE (kerror,*)'Wrong Day of Year value --- >>', doy
          return
      end if

      if(rgma.gt.90.0.or.rgma.lt.-90.0)then
          WRITE (kerror,*)'Wrong GEOMAGNETIC LATITUDE value --- >>', rgma
          return
      end if

C FROM 3-HOURLY TO HOURLY ap (New, interpolates between the three hourly ap values)

      ape(1)=ap(1)
      ape(2)=ap(1)
      ape(38)=ap(13)
      ape(39)=ap(13)

      DO k = 1,13
         ape(k*3-1) = ap(k)
      END DO

      DO k = 1,12
         ape(k*3) = (ap(k)*2 + ap(k+1))/3.0
      END DO

      DO k = 2,13
         ape(k*3-2) = (ap(k-1) + ap(k)*2)/3.0
      END DO

C     TO OBTAIN THE INTEGRAL OF ap.

      k=mod(ut,3)+1

      rap = 0

      DO j = 1,36
         rap = rap + fap(j) * ape(k+j)
      END DO

      if(rap.le.200.)then
         cf=1.0
         return
      end if

      dayno=doy
      if(rgma.lt.0.0)then
         dayno=doy+172
         if(dayno.gt.365)dayno=dayno-365
      end if

      if (dayno.ge.82) then
         rs=(dayno-82.)/45.6+1.
      else
         rs=(dayno+283.)/45.6+1.
      endif
      s1=rs
      facs=rs-s1
      s2=s1+1
      if(s2.eq.9) s2=1

      rgma = abs(rgma)

      rl=(rgma+10.)/20.+1
      if(rl.eq.6.0)rl=5.9
      l1=rl
      facl=rl-l1
      l2=l1+1

C     FACTORS CALCULATIONS

      rap1 = max(300.,rap)
      rap2 = rap1*rap1
      rap3 = rap2*rap1
      n=code(s1,l1)
      cf1=c3(n)*rap3 + c2(n)*rap2 + c1(n)*rap1 + c0(n)
      n=code(s1,l2)
      cf2=c3(n)*rap3 + c2(n)*rap2 + c1(n)*rap1 + c0(n)
      n=code(s2,l1)
      cf3=c3(n)*rap3 + c2(n)*rap2 + c1(n)*rap1 + c0(n)
      n=code(s2,l2)
      cf4=c3(n)*rap3 + c2(n)*rap2 + c1(n)*rap1 + c0(n)

C     INTERPOLATION

      cf = cf1 * (1 - facs) * (1 - facl) + cf2 * (1 - facs) * (facl) +
     &   cf3 * facs * (1 - facl) + cf4 * facs * (facl)

      if (rap.lt.300.) cf = (cf-1.0)*rap/100.-2.*cf+3.

      RETURN
      END
