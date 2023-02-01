c igrf.f, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
C
C Subroutines to compute IGRF parameters for IRI and all functions and
C subroutines required for this computation, including:
C 	igrf_sub, FELDG, FELDCOF, GETSHC, INTERSHC,
C 	EXTRASHC, GEODIP, SPHCAR, GEOMAG, and RECALC
C
C Corrections:
C  1/27/92 Adopted to IGRF-91 coeffcients model
C  2/05/92 Reduce variable names: INTER(P)SHC,EXTRA(P)SHC,INITI(ALI)ZE
C  8/08/95 Updated to IGRF-45-95; new coeff. DGRF90, IGRF95, IGRF95S
C  5/31/00 Updated to IGRF-45-00; new coeff.: IGRF00, IGRF00s
C-Version-mm/dd/yy-Description (Person reporting the correction)
C 2000.01 05/07/01 initial version
C 2000.02 07/11/01 replace feldi(xi,h) by feldi (P. Wilkinson)
C 2000.02 07/11/01 variables EGNR, AGNR,OGNR not used (P. Wilkinson)
c 2000.01 10/28/02 replace TAB/6 blanks, enforce 72/line (D. Simpson)
C 2000.02 11/08/02 change unit for coefficients to 14
C 2000.03 06/05/03 correct DIPL computation (V. Truhlik)
C 2005.00 04/25/05 CALL FELDI and DO 1111 I=1,7 (Alexey Petrov)
C 2005.01 11/10/05 added igrf_dip and geodip (MLAT)
C 2005.02 11/10/05 updated to IGRF-10 version
C 2005.03 12/21/06 GH2(120) -> GH2(144)
C 2007.00 05/18/07 Release of IRI-2007
C Remko Scharroo: call FELDCOF only when integer year changes.
C
      SUBROUTINE igrf_dip(xlat,xlong,year,height,dip,dipl,ymodip)
c-----------------------------------------------------------------------
c INPUT:
c    xlat      geodatic latitude in degrees
c    xlong     geodatic longitude in degrees
c    year      decimal year (year+month/12.0-0.5 or
c                  year+day-of-year/365 or ../366 if leap year)
c    height    height in km
c OUTPUT:
c    dip       magnetic inclination (dip) in degrees
c    dipl      dip latitude in degrees
c    ymodip    modified dip latitude = asin{dip/sqrt[dip^2+cos(LATI)]}
c-----------------------------------------------------------------------

      COMMON /CONST/UMR
      INTEGER IYEAR /0/
      SAVE IYEAR
c
C----------------CALCULATE PROFILES-----------------------------------
c
! Here a little cludge: (by Remko Scharroo)
! Update the coefficients (using FELDCOF) only when the (integer) year
! changes. Saves time and, besides, the magnetic field varies slowly anyhow.
! Also, there was a bug in IRISUB that resulted in FELDCOF to be fed integer
! years before I fixed it, so the outcome is the same with advantage of
! faster processing
      IF (INT(YEAR).NE.IYEAR) THEN
	 IYEAR=INT(YEAR)
      	 CALL FELDCOF(REAL(IYEAR),DIMO)
      ENDIF
      CALL FELDG(XLAT,XLONG,HEIGHT,BNORTH,BEAST,BDOWN,BABS)
      DIP=ASIN(BDOWN/BABS)
      dipdiv=DIP/SQRT(DIP*DIP+cos(XLAT*UMR))
      IF(ABS(dipdiv).GT.1.) dipdiv=SIGN(1.,dipdiv)
      SMODIP=ASIN(dipdiv)
      DIPL=ATAN(BDOWN/2.0/sqrt(BNORTH*BNORTH+BEAST*BEAST))/umr
      YMODIP=SMODIP/UMR
      DIP=DIP/UMR
      RETURN
      END
c
c
C SHELLIG.FOR, Version 2.0, January 1992
C
C  1/27/92 Adopted to IGRF-91 coeffcients model
C  2/05/92 Reduce variable-names: INTER(P)SHC,EXTRA(P)SHC,INITI(ALI)ZE
C  8/08/95 Updated to IGRF-45-95; new coeff. DGRF90, IGRF95, IGRF95S
C  5/31/00 Updated to IGRF-45-00; new coeff.: IGRF00, IGRF00s
C*********************************************************************
C  SUBROUTINES FELDG, FELDCOF, GETSHC, INTERSHC, EXTRASHC
C*********************************************************************
C*********************************************************************
C
C
      SUBROUTINE FELDG(GLAT,GLON,ALT,BNORTH,BEAST,BDOWN,BABS)
c-----------------------------------------------------------------------
C CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
C REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61,
C      1970.
c-----------------------------------------------------------------------
C CHANGES (D. BILITZA, NOV 87):
C   - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA
C   - CALCULATES DIPOL MOMENT
c-----------------------------------------------------------------------
C  INPUT:  ENTRY POINT FELDG
C               GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
C               GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
C               ALT   ALTITUDE IN KM ABOVE SEA LEVEL
C
C          ENTRY POINT FELDC
C               V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
C                       X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
C                       Y-AXIS POINTING TO EQUATOR AT 90 LONG.
C                       Z-AXIS POINTING TO NORTH POLE
C
C          COMMON /MODEL/ AND /IGRF1/ AND /CONST/
C               NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS
C               G(M)    NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
C                       M=NMAX*(NMAX+2)
C               ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN
C                       COORDINATES (6371.2 KM)
C               AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS OF
C                       EARTH ELLIPSOID AS RECOMMENDED BY INTERNAT.
C                       ASTRONOMICAL UNION (6378.160, 6356.775 KM).
C               UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
c-----------------------------------------------------------------------
C  OUTPUT: BABS   MAGNETIC FIELD STRENGTH IN GAUSS
C          BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
C                 TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
C                 POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
C                 AND DOWNWARD.
C-----------------------------------------------------------------------
      DIMENSION         B(3)
      COMMON/IGRF/      XI(3),H(196)
      COMMON/MODEL/     NMAX,G(196)
      COMMON/IGRF1/     ERA,AQUAD,BQUAD
      COMMON/CONST/	UMR
C
C-- IS RECORDS ENTRY POINT
C
C*****ENTRY POINT  FELDG  TO BE USED WITH GEODETIC CO-ORDINATES
      IS=1
      RLAT=GLAT*UMR
      CT=SIN(RLAT)
      ST=COS(RLAT)
      D=SQRT(AQUAD-(AQUAD-BQUAD)*CT*CT)
      RLON=GLON*UMR
      CP=COS(RLON)
      SP=SIN(RLON)
      ZZZ=(ALT+BQUAD/D)*CT/ERA
      RHO=(ALT+AQUAD/D)*ST/ERA
      XXX=RHO*CP
      YYY=RHO*SP
      RQ=1./(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      XI(1)=XXX*RQ
      XI(2)=YYY*RQ
      XI(3)=ZZZ*RQ
      IHMAX=NMAX*NMAX+1
      LAST=IHMAX+NMAX+NMAX
      IMAX=NMAX+NMAX-1
      DO 8 I=IHMAX,LAST
8     H(I)=G(I)
      DO 6 K=1,3,2
      I=IMAX
      IH=IHMAX
1     IL=IH-I
      F=2./FLOAT(I-K+2)
      X=XI(1)*F
      Y=XI(2)*F
      Z=XI(3)*(F+F)
      I=I-2
      IF(I-1)5,4,2
2     DO 3 M=3,I,2
      H(IL+M+1)=G(IL+M+1)+Z*H(IH+M+1)+X*(H(IH+M+3)-H(IH+M-1))
     A                               -Y*(H(IH+M+2)+H(IH+M-2))
3     H(IL+M)=G(IL+M)+Z*H(IH+M)+X*(H(IH+M+2)-H(IH+M-2))
     A                         +Y*(H(IH+M+3)+H(IH+M-1))
4     H(IL+2)=G(IL+2)+Z*H(IH+2)+X*H(IH+4)-Y*(H(IH+3)+H(IH))
      H(IL+1)=G(IL+1)+Z*H(IH+1)+Y*H(IH+4)+X*(H(IH+3)-H(IH))
5     H(IL)=G(IL)+Z*H(IH)+2.*(X*H(IH+1)+Y*H(IH+2))
      IH=IL
      IF(I.GE.K)GOTO 1
6     CONTINUE
      IF(IS.EQ.3)RETURN
      S=.5*H(1)+2.*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      T=(RQ+RQ)*SQRT(RQ)
      BXXX=T*(H(3)-S*XXX)
      BYYY=T*(H(4)-S*YYY)
      BZZZ=T*(H(2)-S*ZZZ)
      IF(IS.EQ.2)GOTO 7
      BABS=SQRT(BXXX*BXXX+BYYY*BYYY+BZZZ*BZZZ)
      BEAST=BYYY*CP-BXXX*SP
      BRHO=BYYY*SP+BXXX*CP
      BNORTH=BZZZ*ST-BRHO*CT
      BDOWN=-BZZZ*CT-BRHO*ST
      RETURN
7     B(1)=BXXX
      B(2)=BYYY
      B(3)=BZZZ
      RETURN
      END
C
C
        SUBROUTINE FELDCOF(YEAR,DIMO)
c-----------------------------------------------------------------------
C  DETERMINES COEFFICIENTS AND DIPOL MOMENT FROM IGRF MODELS
C
C       INPUT:  YEAR    DECIMAL YEAR FOR WHICH GEOMAGNETIC FIELD IS TO
C                       BE CALCULATED
C       OUTPUT: DIMO    GEOMAGNETIC DIPOL MOMENT IN GAUSS (NORMALIZED
C                       TO EARTH'S RADIUS) AT THE TIME (YEAR)
C  D. BILITZA, NSSDC, GSFC, CODE 633, GREENBELT, MD 20771,
C       (301)286-9536   NOV 1987.
C  ### updated to IGRF-10 version -dkb- 11/10/2005
c-----------------------------------------------------------------------
      PARAMETER	(NUMYE=14)
C ### NUMYE = number of years represented by IGRF + one with rate
      CHARACTER*13    FILMOD(NUMYE+1)
      DIMENSION       GH1(196),GH2(196),GHA(196),GHB(196)
      DOUBLE PRECISION X,F0,F
      COMMON/MODEL/   NMAX,GHB
      COMMON/IGRF1/   ERA,AQUAD,BQUAD
      COMMON/CONST/   UMR
C ### updated to 2010
        DATA  FILMOD   / 'dgrf1945.dat','dgrf1950.dat','dgrf1955.dat',           
     1    'dgrf1960.dat','dgrf1965.dat','dgrf1970.dat','dgrf1975.dat',
     2    'dgrf1980.dat','dgrf1985.dat','dgrf1990.dat','dgrf1995.dat',
     3    'dgrf2000.dat','dgrf2005.dat','igrf2010.dat','igrf2010s.dat'/
      DATA ERA/6371.2/, AQUAD/40680924.9856/, BQUAD/40408588.4006/
      SAVE
C
C
C  IS=0 FOR SCHMIDT NORMALIZATION   IS=1 GAUSS NORMALIZATION
C
        IS = 0
	IER = 0
C-- DETERMINE IGRF-YEARS FOR INPUT-YEAR
        L = (INT(YEAR) - 1945)/5 + 1
        IF (L.LT.1) L=1
        IF (L.GT.NUMYE) L=NUMYE
        DTE1 = L * 5. + 1940.
        DTE2 = DTE1 + 5.
C-- GET IGRF COEFFICIENTS FOR THE BOUNDARY YEARS
	CALL GETSHC (FILMOD(L), NMAX1, ERA, GH1, IER)
        IF (IER .NE. 0) STOP "IRI2007: Error reading IGRF file"
        CALL GETSHC (FILMOD(L+1), NMAX2, ERA, GH2, IER)
        IF (IER .NE. 0) STOP "IRI2007: Error reading IGRF file"
        IF (L.LT.NUMYE) THEN
C-- DETERMINE IGRF COEFFICIENTS FOR YEAR
          CALL INTERSHC (YEAR, DTE1, NMAX1, GH1, DTE2,
     1          NMAX2, GH2, NMAX, GHA)
        ELSE
          CALL EXTRASHC (YEAR, DTE1, NMAX1, GH1, NMAX2,
     1          GH2, NMAX, GHA)
        ENDIF
C-- DETERMINE MAGNETIC DIPOL MOMENT AND COEFFICIENTS G
        F0=0.D0
        DO 1234 J=1,3
           F = GHA(J) * 1.D-5
           F0 = F0 + F * F
1234    CONTINUE
        DIMO = DSQRT(F0)

        GHB(1) =  0.0
        I=2
        F0=1.D-5
        IF(IS.EQ.0) F0=-F0
        SQRT2=SQRT(2.)

      DO 9 N=1,NMAX
        X = N
        F0 = F0 * X * X / (4.D0 * X - 2.D0)
        IF(IS.EQ.0) F0 = F0 * (2.D0 * X - 1.D0) / X
        F = F0 * 0.5D0
        IF(IS.EQ.0) F = F * SQRT2
        GHB(I) = GHA(I-1) * F0
        I = I+1
      DO 9 M=1,N
        F = F * (X + M) / (X - M + 1.D0)
        IF(IS.EQ.0) F = F * DSQRT((X - M + 1.D0) / (X + M))
        GHB(I) = GHA(I-1) * F
        GHB(I+1) = GHA(I) * F
        I=I+2
9     CONTINUE
        RETURN
        END
C
C
        SUBROUTINE GETSHC (FSPEC, NMAX, ERA, GH, IER)

C ===============================================================
C
C       Version 1.01
C
C       Reads spherical harmonic coefficients from the specified
C       file into an array.
C
C       Input:
C           FSPEC - File specification
C
C       Output:
C           NMAX  - Maximum degree and order of model
C           ERA   - Earth's radius associated with the spherical
C                   harmonic coefficients, in the same units as
C                   elevation
C           GH    - Schmidt quasi-normal internal spherical
C                   harmonic coefficients
C           IER   - Error number: =  0, no error
C                                 = -2, records out of order
C                                 = FORTRAN run-time error number
C
C       A. Zunde
C       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
C
C ===============================================================

      CHARACTER  FSPEC*(*), dirnam*80
      integer	freeunit
      DIMENSION       GH(*)
      common	/iounit/konsol,kerror
        
        do 1 j=1,196  
1          GH(j)=0.0
C ---------------------------------------------------------------
C       Open coefficient file. Read past first header record.
C       Read degree and order of model and Earth radius.
C ---------------------------------------------------------------
      IU=freeunit()
      dirnam='/user/altim'
      call checkenv('ALTIM',dirnam,l)
      IER=0
      dirnam(l+1:)='/data/iri2007/'//fspec
      if (konsol.ge.0) write (konsol,*) "Reading: ",dirnam
      OPEN (IU, FILE=dirnam,STATUS='OLD', IOSTAT=IER, ERR=999)

      READ (IU, *, IOSTAT=IER, ERR=999)
      READ (IU, *, IOSTAT=IER, ERR=999) NMAX, ERA, XMYEAR
        nm=nmax*(nmax+2)                
        READ (IU, *, IOSTAT=IER, ERR=999) (GH(i),i=1,nm) 
999     CLOSE (IU)
        RETURN
        END
C
C
        SUBROUTINE INTERSHC (DATE, DTE1, NMAX1, GH1, DTE2,
     1                        NMAX2, GH2, NMAX, GH)

C ===============================================================
C
C       Version 1.01
C
C       Interpolates linearly, in time, between two spherical
C       harmonic models.
C
C       Input:
C           DATE  - Date of resulting model (in decimal year)
C           DTE1  - Date of earlier model
C           NMAX1 - Maximum degree and order of earlier model
C           GH1   - Schmidt quasi-normal internal spherical
C                   harmonic coefficients of earlier model
C           DTE2  - Date of later model
C           NMAX2 - Maximum degree and order of later model
C           GH2   - Schmidt quasi-normal internal spherical
C                   harmonic coefficients of later model
C
C       Output:
C           GH    - Coefficients of resulting model
C           NMAX  - Maximum degree and order of resulting model
C
C       A. Zunde
C       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
C
C ===============================================================

        DIMENSION       GH1(*), GH2(*), GH(*)

C ---------------------------------------------------------------
C       The coefficients (GH) of the resulting model, at date
C       DATE, are computed by linearly interpolating between the
C       coefficients of the earlier model (GH1), at date DTE1,
C       and those of the later model (GH2), at date DTE2. If one
C       model is smaller than the other, the interpolation is
C       performed with the missing coefficients assumed to be 0.
C ---------------------------------------------------------------

        FACTOR = (DATE - DTE1) / (DTE2 - DTE1)

        IF (NMAX1 .EQ. NMAX2) THEN
            K = NMAX1 * (NMAX1 + 2)
            NMAX = NMAX1
        ELSE IF (NMAX1 .GT. NMAX2) THEN
            K = NMAX2 * (NMAX2 + 2)
            L = NMAX1 * (NMAX1 + 2)
            DO 1122 I = K + 1, L
1122            GH(I) = GH1(I) + FACTOR * (-GH1(I))
            NMAX = NMAX1
        ELSE
            K = NMAX1 * (NMAX1 + 2)
            L = NMAX2 * (NMAX2 + 2)
            DO 1133 I = K + 1, L
1133            GH(I) = FACTOR * GH2(I)
            NMAX = NMAX2
        ENDIF

        DO 1144 I = 1, K
1144        GH(I) = GH1(I) + FACTOR * (GH2(I) - GH1(I))

        RETURN
        END
C
C
        SUBROUTINE EXTRASHC (DATE, DTE1, NMAX1, GH1, NMAX2,
     1                        GH2, NMAX, GH)

C ===============================================================
C
C       Version 1.01
C
C       Extrapolates linearly a spherical harmonic model with a
C       rate-of-change model.
C
C       Input:
C           DATE  - Date of resulting model (in decimal year)
C           DTE1  - Date of base model
C           NMAX1 - Maximum degree and order of base model
C           GH1   - Schmidt quasi-normal internal spherical
C                   harmonic coefficients of base model
C           NMAX2 - Maximum degree and order of rate-of-change
C                   model
C           GH2   - Schmidt quasi-normal internal spherical
C                   harmonic coefficients of rate-of-change model
C
C       Output:
C           GH    - Coefficients of resulting model
C           NMAX  - Maximum degree and order of resulting model
C
C       A. Zunde
C       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
C
C ===============================================================

        DIMENSION       GH1(*), GH2(*), GH(*)

C ---------------------------------------------------------------
C       The coefficients (GH) of the resulting model, at date
C       DATE, are computed by linearly extrapolating the coef-
C       ficients of the base model (GH1), at date DTE1, using
C       those of the rate-of-change model (GH2), at date DTE2. If
C       one model is smaller than the other, the extrapolation is
C       performed with the missing coefficients assumed to be 0.
C ---------------------------------------------------------------

        FACTOR = (DATE - DTE1)

        IF (NMAX1 .EQ. NMAX2) THEN
            K = NMAX1 * (NMAX1 + 2)
            NMAX = NMAX1
        ELSE IF (NMAX1 .GT. NMAX2) THEN
            K = NMAX2 * (NMAX2 + 2)
            L = NMAX1 * (NMAX1 + 2)
            DO 1155 I = K + 1, L
1155            GH(I) = GH1(I)
            NMAX = NMAX1
        ELSE
            K = NMAX1 * (NMAX1 + 2)
            L = NMAX2 * (NMAX2 + 2)
            DO 1166 I = K + 1, L
1166            GH(I) = FACTOR * GH2(I)
            NMAX = NMAX2
        ENDIF

        DO 1177 I = 1, K
1177        GH(I) = GH1(I) + FACTOR * GH2(I)

        RETURN
        END
C
C
      SUBROUTINE GEODIP(IYR,SLA,SLO,DLA,DLO,J)

C  Calculates dipole geomagnetic coordinates from geocentric coordinates
C  or vice versa.

C                     J=0           J=1
C		INPUT:     J,SLA,SLO     J,DLA,DLO
C		OUTPUT:     DLA,DLO       SLA,SLO

C  Last revision: November 2005 (Vladimir Papitashvili)
C  The code is modifed from GEOCOR written by V.Popov and V.Papitashvili
C  in mid-1980s.

         COMMON /CONST/UMR

C  Earth radius (km) RE = 6371.2

C  The radius of the sphere to compute the coordinates (in Re)
C        RH = (RE + HI)/RE
         R = 1.

         if(j.gt.0) goto 1234

         COL = (90.- SLA)*UMR
         RLO = SLO*UMR
      CALL SPHCAR(R,COL,RLO,X,Y,Z,1)
      CALL GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
      CALL SPHCAR(RM,TH,PF,XM,YM,ZM,-1)
         SZM = ZM
         DLO = PF/UMR
         DCO = TH/UMR
         DLA = 90.- DCO
         RETURN

1234     continue
      COL = (90.- DLA)*UMR
      RLO = DLO*UMR
      CALL SPHCAR(R,COL,RLO,XM,YM,ZM,1)
      CALL GEOMAG(X,Y,Z,XM,YM,ZM,-1,IYR)
      CALL SPHCAR(RM,TH,PF,X,Y,Z,-1)
        SZM = ZM
        SLO = PF/UMR
        SCO = TH/UMR
        SLA = 90.- SCO

      RETURN
      END

C  *********************************************************************

      SUBROUTINE SPHCAR(R,TETA,PHI,X,Y,Z,J)

C   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
C    (TETA AND PHI IN RADIANS).

C                  J>0            J<0
C-----INPUT:   J,R,TETA,PHI     J,X,Y,Z
C----OUTPUT:      X,Y,Z        R,TETA,PHI

C  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
C      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
C      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)

        IMPLICIT NONE

        REAL R,TETA,PHI,X,Y,Z,SQ

        INTEGER J

      IF(J.GT.0) GOTO 3
      SQ=X**2+Y**2
      R=SQRT(SQ+Z**2)
      IF (SQ.NE.0.) GOTO 2
      PHI=0.
      IF (Z.LT.0.) GOTO 1
      TETA=0.
      RETURN
  1   TETA=3.141592654
      RETURN
  2   SQ=SQRT(SQ)
      PHI=ATAN2(Y,X)
      TETA=ATAN2(SQ,Z)
      IF (PHI.LT.0.) PHI=PHI+6.28318531
      RETURN
  3   SQ=R*SIN(TETA)
      X=SQ*COS(PHI)
      Y=SQ*SIN(PHI)
      Z=R*COS(TETA)

      RETURN
      END

C  *********************************************************************

      SUBROUTINE GEOMAG(XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J,IYR)

C CONVERTS GEOCENTRIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.
C IYR IS YEAR NUMBER (FOUR DIGITS).

C                           J>0                J<0
C-----INPUT:  J,XGEO,YGEO,ZGEO,IYR   J,XMAG,YMAG,ZMAG,IYR
C-----OUTPUT:    XMAG,YMAG,ZMAG        XGEO,YGEO,ZGEO

C  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
C      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
C      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)

      IMPLICIT NONE

      REAL XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,ST0,CT0,SL0,CL0,CTCL,
     *       STCL,CTSL,STSL

      INTEGER J,IYR,II

      COMMON/C1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL
      DATA II/1/
      SAVE II
      IF (IYR.NE.II) THEN
         II=IYR
         CALL RECALC(II,0)
      ENDIF
      IF (J.LT.0) THEN
         XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
         YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
         ZGEO=ZMAG*CT0-XMAG*ST0
      ELSE
         XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
         YMAG=YGEO*CL0-XGEO*SL0
         ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
      ENDIF

      RETURN
      END

C  *********************************************************************

      SUBROUTINE RECALC(IYR,IDAY)

C  Modified to accept years from 1900 through 2015 using the DGRF &
C     IGRF-10 coeficients (updated by V. Papitashvili, November 2005)
C     IYR = YEAR NUMBER (FOUR DIGITS)
C     IDAY = DAY NUMBER

      IMPLICIT NONE

      REAL ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,
     2       F2,F1,G10,G11,H11,SQQ,SQR

      INTEGER IYR,IDAY

      COMMON/C1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL

C  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
C  VALUES FOR THE NEAREST EPOCHS:

       IF (IYR.LT.1905) THEN                             !1900-1905
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1900.)/5.
           F1=1.D0-F2
           G10=31543.*F1+31464.*F2
           G11=-2298.*F1-2298.*F2
           H11= 5922.*F1+5909.*F2
       ELSEIF (IYR.LT.1910) THEN                         !1905-1910
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1905.)/5.
           F1=1.D0-F2
           G10=31464.*F1+31354.*F2
           G11=-2298.*F1-2297.*F2
           H11= 5909.*F1+5898.*F2
       ELSEIF (IYR.LT.1915) THEN                         !1910-1915
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1910.)/5.
           F1=1.D0-F2
           G10=31354.*F1+31212.*F2
           G11=-2297.*F1-2306.*F2
           H11= 5898.*F1+5875.*F2
       ELSEIF (IYR.LT.1920) THEN                         !1915-1920
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1915.)/5.
           F1=1.D0-F2
           G10=31212.*F1+31060.*F2
           G11=-2306.*F1-2317.*F2
           H11= 5875.*F1+5845.*F2
       ELSEIF (IYR.LT.1925) THEN                         !1920-1925
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1920.)/5.
           F1=1.D0-F2
           G10=31060.*F1+30926.*F2
           G11=-2317.*F1-2318.*F2
           H11= 5845.*F1+5817.*F2
       ELSEIF (IYR.LT.1930) THEN                         !1925-1930
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1925.)/5.
           F1=1.D0-F2
           G10=30926.*F1+30805.*F2
           G11=-2318.*F1-2316.*F2
           H11= 5817.*F1+5808.*F2
        ELSEIF (IYR.LT.1935) THEN                        !1930-1935
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1930.)/5.
           F1=1.D0-F2
           G10=30805.*F1+30715.*F2
           G11=-2316.*F1-2306.*F2
           H11= 5808.*F1+5812.*F2
        ELSEIF (IYR.LT.1940) THEN                        !1935-1940
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1935.)/5.
           F1=1.D0-F2
           G10=30715.*F1+30654.*F2
           G11=-2306.*F1-2292.*F2
           H11= 5812.*F1+5821.*F2
        ELSEIF (IYR.LT.1945) THEN                        !1940-1945
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1940.)/5.
           F1=1.D0-F2
           G10=30654.*F1+30594.*F2
           G11=-2292.*F1-2285.*F2
           H11= 5821.*F1+5810.*F2
        ELSEIF (IYR.LT.1950) THEN                        !1945-1950
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1945.)/5.
           F1=1.D0-F2
           G10=30594.*F1+30554.*F2
           G11=-2285.*F1-2250.*F2
           H11= 5810.*F1+5815.*F2
        ELSEIF (IYR.LT.1955) THEN                        !1950-1955
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1950.)/5.
           F1=1.D0-F2
           G10=30554.*F1+30500.*F2
           G11=-2250.*F1-2215.*F2
           H11= 5815.*F1+5820.*F2
        ELSEIF (IYR.LT.1960) THEN                        !1955-1960
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1955.)/5.
           F1=1.D0-F2
           G10=30500.*F1+30421.*F2
           G11=-2215.*F1-2169.*F2
           H11= 5820.*F1+5791.*F2
        ELSEIF (IYR.LT.1965) THEN                        !1960-1965
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1960.)/5.
           F1=1.D0-F2
           G10=30421.*F1+30334.*F2
           G11=-2169.*F1-2119.*F2
           H11= 5791.*F1+5776.*F2
        ELSEIF (IYR.LT.1970) THEN                        !1965-1970
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1965.)/5.
           F1=1.D0-F2
           G10=30334.*F1+30220.*F2
           G11=-2119.*F1-2068.*F2
           H11= 5776.*F1+5737.*F2
        ELSEIF (IYR.LT.1975) THEN                        !1970-1975
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1970.)/5.
           F1=1.D0-F2
           G10=30220.*F1+30100.*F2
           G11=-2068.*F1-2013.*F2
           H11= 5737.*F1+5675.*F2
        ELSEIF (IYR.LT.1980) THEN                        !1975-1980
           F2=(DFLOAT(IYR)+DFLOAT(IDAY)/365.-1975.)/5.
           F1=1.D0-F2
           G10=30100.*F1+29992.*F2
           G11=-2013.*F1-1956.*F2
           H11= 5675.*F1+5604.*F2
        ELSEIF (IYR.LT.1985) THEN                        !1980-1985
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1980.)/5.
           F1=1.D0-F2
           G10=29992.*F1+29873.*F2
           G11=-1956.*F1-1905.*F2
           H11= 5604.*F1+5500.*F2
        ELSEIF (IYR.LT.1990) THEN                        !1985-1990
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1985.)/5.
           F1=1.D0-F2
           G10=29873.*F1+29775.*F2
           G11=-1905.*F1-1848.*F2
           H11= 5500.*F1+5406.*F2
        ELSEIF (IYR.LT.1995) THEN                        !1990-1995
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1990.)/5.
           F1=1.D0-F2
           G10=29775.*F1+29692.*F2
           G11=-1848.*F1-1784.*F2
           H11= 5406.*F1+5306.*F2
        ELSEIF (IYR.LT.2000) THEN                        !1995-2000
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-1995.)/5.
           F1=1.D0-F2
           G10=29692.*F1+29619.4*F2
           G11=-1784.*F1-1728.2*F2
           H11= 5306.*F1+5186.1*F2
        ELSEIF (IYR.LT.2005) THEN                        !2000-2005
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-2000.)/5.
           F1=1.D0-F2
           G10=29619.4*F1+29554.63*F2
           G11=-1728.2*F1-1669.05*F2
           H11= 5186.1*F1+5077.99*F2
        ELSEIF (IYR.LT.2010) THEN                        !2005-2010
           F2=(FLOAT(IYR)+FLOAT(IDAY)/365.-2005.)/5.
           F1=1.D0-F2
           G10=29554.63*F1+29496.5*F2
           G11=-1669.05*F1-1585.9*F2
           H11= 5077.99*F1+4945.1*F2
        ELSE                                            !2010-2015
           F2=FLOAT(IYR)+FLOAT(IDAY)/365.-2010.
           G10=29496.5-11.4*F2
           G11=-1585.9+16.7*F2
           H11= 4945.1-28.8*F2
        ENDIF

C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD
C  SYSTEM:
C  SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
C         ST0 * CL0                ST0 * SL0                CT0

      SQQ=G11**2+H11**2
      SQR=SQRT(G10**2+SQQ)
      SQQ=SQRT(SQQ)
      SL0=-H11/SQQ
      CL0=-G11/SQQ
      ST0=SQQ/SQR
      CT0=G10/SQR
      STCL=ST0*CL0
      STSL=ST0*SL0
      CTSL=CT0*SL0
      CTCL=CT0*CL0

C  THE CALCULATIONS ARE TERMINATED IF ONLY GEO-MAG TRANSFORMATION
C  IS TO BE DONE

      RETURN
      END
