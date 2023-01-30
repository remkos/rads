      SUBROUTINE PERTH2( DLAT,DLON,TIME,TIDE,ISDATA )
*
*  Name - PERTH                  Name derivation - PREdict Tidal Heights
*
*  Function -  to compute the ocean tidal height at a given time
*              and location from grids of harmonic constants.
*              Current version uses the 8 largest constituents in the
*              semidiurnal & diurnal bands, with other tides inferred.
*              (Long period tides are NOT computed by this routine.)
*
*  Language - Fortran 77
*
*  Arguments -
*     name      type  I/O               description
*     ----      ----  ---               -----------
*     DLAT       D     I    north latitude (in degrees) for desired
*                           location.
*
*     DLON       D     I    east longitude (in degrees).
*
*     TIME       D     I    desired UTC time, in (decimal) Modified Julian Date
*                           e.g., 1 Jan 1990 noon = 47892.5
*
*     TIDE       D     O    computed tidal height.  The units will be
*                           the same as the 'amplitude' units on the
*                           input tidal grids (usually cm).
*
*     ISDATA     L     O    logical denoting whether tide data exist at
*                           desired location.  If FALSE, then TIDE is
*                           not modified.
*
*
*  Usage notes -
*     All 8 input tide files should be concatenated, in order of frequency,
*     e.g., do:
*         cat q1.d o1.d p1.d k1.d n2.d m2.d s2.d k2.d > fort.30
*     (or the equivalent for your operating system).
*
*  Processing logic -
*     Tidal constants at desired location are found by bilinear
*     interpolation from input grid files.  The astronomical mean
*     longitudes (determined by ASTROL) are most accurate for the
*     period 1990 - 2010.  Nodal corrections are applied to all
*     lunar tides.  Sixteen minor tides are inferred from major tides.
*
*  File references -
*     Input tidal grids are read on first call via unit LU (set in DATA).
*
*  Important local variables -
*     ARRAY - holds inphase & quadrature tidal constants as follows:
*        ARRAY(i,j,k,l) where
*        i,j are longitude/latitude indices.
*        k = 1,2  for Hcos(G) or Hsin(G).
*        l = 1,...,8 (max) - for each tidal constituent.
*     SHPN - holds the mean astronomical longitudes of interest;
*        this array is equivalenced to S, H, P, & OMEGA
*
*     Programming note: This routine is written for general use, which
*     is not necessarily efficient for all applications.
*     If it is desired to compute tidal heights at, e.g., every 1-sec
*     observation along an arc, the program could be speeded up by
*     use of several approximations.  Contact author for details.
*     For machines like the Cray, the DOUBLE PRECISION should be removed.
*
*  Error processing -
*     An attempt is made is ensure that the input grids are read in the
*     correct order; this procedure looks for the constituent names
*     in the TITLE part of the input files.
*
*  Technical references -
*     A. T. Doodson & H. Warburg, Admiralty Manual of Tides, HMSO, 1941.
*
*  History -
*   version   date    programmer        change description
*   -------   ----    ----------        ------------------
*     1.0    6/03/93     R Ray    Initial version.
*     1.1    1/10/95     R Ray    Fixed missing declaration of NINETY.
*     1.2    3/15/95     R Ray    TLOAD ensures consistent nulls in all tides.
*     1.3    5/05/95     R Ray    Keep last ISDATA internally.
*     2.0    5/19/95     R Ray    Infer minor tides; remove SELECT argument.
*                                 Bypass nodal calculations if time < 30 days.
*     2.1    1/13/98     R Ray    Fixed bad argument to "pi1" tide.
*
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  NINETY
      PARAMETER (MX=720, MY=361, NG=8, NT=26)
      LOGICAL    INIT,ISDATA,SELECT(NG),WRAP,ISDATA0
      REAL       ARRAY(MX,MY,2,NG)
      REAL       LATMIN,LATMAX,LONMIN,LONMAX,UNDEF,H12(2,NT)
      REAL       SINN,SIN2N,COSN,COS2N
      DIMENSION  SHPN(4), ARG(NT), F(NT), U(NT)
      EQUIVALENCE (SHPN(1),S),(SHPN(2),H),(SHPN(3),P),(SHPN(4),OMEGA)
      SAVE
      PARAMETER (UNDEF=99999.0, PI=3.141592654D0, RAD=PI/180.D0)
      PARAMETER (TWO=2.D0, THREE=3.D0, FOUR=4.D0)
      PARAMETER (FIFTEN=15.D0, THIRTY=30.D0, NINETY=90.D0, PP=282.8D0)
      DATA       DLAT0,DLON0,TIME0,TIME1/4*999.d0/
      DATA       INIT/.TRUE./, LU/30/, SELECT/8*.TRUE./

*     on first call, read grids of harmonic constants
*     -----------------------------------------------
      IF (INIT) THEN
         INIT = .FALSE.
         CALL TLOAD( ARRAY,MX,MY, SELECT, LU, UNDEF,
     *               NX,NY,LATMIN,LATMAX,LONMIN,LONMAX,NGRIDS )
         DX = (LONMAX - LONMIN)/(NX - 1)
         WRAP = ABS(LONMAX - LONMIN - 360.) .LT. 2.0*DX
         DO 10 I=1,NT
            F(I) = 1.D0
            U(I) = 0.D0
 10      CONTINUE
      ENDIF
*
*     determine tidal constants at this location
*     ------------------------------------------
      IF (DLAT.NE.DLAT0 .OR. DLON.NE.DLON0) THEN
          CALL GRSINT( ARRAY,MX,MY,2*NGRIDS,NX,NY,UNDEF,
     *                 LATMIN,LATMAX,LONMIN,LONMAX,WRAP,
     *                 REAL(DLAT),REAL(DLON),H12,ISDATA )
          DLAT0 = DLAT
          DLON0 = DLON
          ISDATA0 = ISDATA
*
*         infer minor tides at this location
*         ----------------------------------
          H12(1, 9) = 0.263 *H12(1,1) - 0.0252*H12(1,2)
          H12(2, 9) = 0.263 *H12(2,1) - 0.0252*H12(2,2)
          H12(1,10) = 0.297 *H12(1,1) - 0.0264*H12(1,2)
          H12(2,10) = 0.297 *H12(2,1) - 0.0264*H12(2,2)
          H12(1,11) = 0.164 *H12(1,1) + 0.0048*H12(1,2)
          H12(2,11) = 0.164 *H12(2,1) + 0.0048*H12(2,2)
          H12(1,12) = 0.0140*H12(1,2) + 0.0101*H12(1,4)
          H12(2,12) = 0.0140*H12(2,2) + 0.0101*H12(2,4)
          H12(1,13) = 0.0389*H12(1,2) + 0.0282*H12(1,4)
          H12(2,13) = 0.0389*H12(2,2) + 0.0282*H12(2,4)
          H12(1,14) = 0.0064*H12(1,2) + 0.0060*H12(1,4)
          H12(2,14) = 0.0064*H12(2,2) + 0.0060*H12(2,4)
          H12(1,15) = 0.0030*H12(1,2) + 0.0171*H12(1,4)
          H12(2,15) = 0.0030*H12(2,2) + 0.0171*H12(2,4)
          H12(1,16) =-0.0015*H12(1,2) + 0.0152*H12(1,4)
          H12(2,16) =-0.0015*H12(2,2) + 0.0152*H12(2,4)
          H12(1,17) =-0.0065*H12(1,2) + 0.0155*H12(1,4)
          H12(2,17) =-0.0065*H12(2,2) + 0.0155*H12(2,4)
          H12(1,18) =-0.0389*H12(1,2) + 0.0836*H12(1,4)
          H12(2,18) =-0.0389*H12(2,2) + 0.0836*H12(2,4)
          H12(1,19) =-0.0431*H12(1,2) + 0.0613*H12(1,4)
          H12(2,19) =-0.0431*H12(2,2) + 0.0613*H12(2,4)
          H12(1,20) = 0.264 *H12(1,5) - 0.0253*H12(1,6)
          H12(2,20) = 0.264 *H12(2,5) - 0.0253*H12(2,6)
          H12(1,21) = 0.298 *H12(1,5) - 0.0264*H12(1,6)
          H12(2,21) = 0.298 *H12(2,5) - 0.0264*H12(2,6)
          H12(1,22) = 0.165 *H12(1,5) + 0.00487*H12(1,6)
          H12(2,22) = 0.165 *H12(2,5) + 0.00487*H12(2,6)
          H12(1,23) = 0.0040*H12(1,6) + 0.0074*H12(1,7)
          H12(2,23) = 0.0040*H12(2,6) + 0.0074*H12(2,7)
          H12(1,24) = 0.0131*H12(1,6) + 0.0326*H12(1,7)
          H12(2,24) = 0.0131*H12(2,6) + 0.0326*H12(2,7)
          H12(1,25) = 0.0033*H12(1,6) + 0.0082*H12(1,7)
          H12(2,25) = 0.0033*H12(2,6) + 0.0082*H12(2,7)
          H12(1,26) = 0.0585*H12(1,7)
          H12(2,26) = 0.0585*H12(2,7)
      ELSE
          ISDATA = ISDATA0
      ENDIF
      IF (.NOT.ISDATA) RETURN
*
*     determine equilibrium tidal arguments
*     -------------------------------------
      IF (TIME.NE.TIME0) THEN
         TIME0 = TIME
         HOUR = (TIME - INT(TIME))*24.D0
         T1 = FIFTEN*HOUR
         T2 = THIRTY*HOUR
         CALL ASTROL( TIME, SHPN )

         ARG(1) = T1 + H - THREE*S + P - NINETY   ! Q1
         ARG(2) = T1 + H - TWO*S - NINETY         ! O1
         ARG(3) = T1 - H - NINETY                 ! P1
         ARG(4) = T1 + H + NINETY                 ! K1
         ARG(5) = T2 + TWO*H - THREE*S + P        ! N2
         ARG(6) = T2 + TWO*H - TWO*S              ! M2
         ARG(7) = T2                              ! S2
         ARG(8) = T2 + TWO*H                      ! K2
         ARG( 9) = T1 - FOUR*S + H + TWO*P - NINETY    ! 2Q1
         ARG(10) = T1 - FOUR*S + THREE*H - NINETY      ! sigma1
         ARG(11) = T1 - THREE*S + THREE*H - P - NINETY ! rho1
         ARG(12) = T1 - S + H - P + NINETY             ! M1
         ARG(13) = T1 - S + H + P + NINETY             ! M1
         ARG(14) = T1 - S + THREE*H - P + NINETY       ! chi1
         ARG(15) = T1 - TWO*H + PP - NINETY            ! pi1
         ARG(16) = T1 + THREE*H + NINETY               ! phi1
         ARG(17) = T1 + S - H + P + NINETY             ! theta1
         ARG(18) = T1 + S + H - P + NINETY             ! J1
         ARG(19) = T1 + TWO*S + H + NINETY             ! OO1
         ARG(20) = T2 - FOUR*S + TWO*H + TWO*P         ! 2N2
         ARG(21) = T2 - FOUR*S + FOUR*H                ! mu2
         ARG(22) = T2 - THREE*S + FOUR*H - P           ! nu2
         ARG(23) = T2 - S + P + 180.D0                 ! lambda2
         ARG(24) = T2 - S + TWO*H - P + 180.D0         ! L2
         ARG(25) = T2 - S + TWO*H + P                  ! L2
         ARG(26) = T2 - H + PP                         ! T2
      ENDIF
*
*     determine nodal corrections f and u    
*        Note: Update this code next iteration of model  -RDR
*     -----------------------------------
      IF (ABS(TIME-TIME1).GT.30.D0) THEN
         TIME1 = TIME
         SINN = SIN(OMEGA*RAD)
         COSN = COS(OMEGA*RAD)
         SIN2N = SIN(TWO*OMEGA*RAD)
         COS2N = COS(TWO*OMEGA*RAD)

         F(1) = 1.009 + 0.187*COSN - 0.015*COS2N
         F(2) = F(1)
         F(4) = 1.006 + 0.115*COSN - 0.009*COS2N
         F(5) = 1.000 - 0.037*COSN
         F(6) = F(5)
         F(8) = 1.024 + 0.286*COSN + 0.008*COS2N
         F( 9) = SQRT((1.0 + 0.189*COSN - 0.0058*COS2N)**2 +
     *                (0.189*SINN - 0.0058*SIN2N)**2)
         F(10) = F(9)
         F(11) = F(9)
         F(12) = SQRT((1.0 + 0.185*COSN)**2 + (0.185*SINN)**2)
         F(13) = SQRT((1.0 + 0.201*COSN)**2 + (0.201*SINN)**2)
         F(14) = SQRT((1.0 + 0.221*COSN)**2 + (0.221*SINN)**2)
         F(18) = SQRT((1.0 + 0.198*COSN)**2 + (0.198*SINN)**2)
         F(19) = SQRT((1.0 + 0.640*COSN + 0.134*COS2N)**2 +
     *                (0.640*SINN + 0.134*SIN2N)**2 )
         F(20) = SQRT((1.0 - 0.0373*COSN)**2 + (0.0373*SINN)**2)
         F(21) = F(20)
         F(22) = F(20)
         F(24) = F(20)
         F(25) = SQRT((1.0 + 0.441*COSN)**2 + (0.441*SINN)**2)

         U(1) = 10.8*SINN - 1.3*SIN2N
         U(2) = U(1)
         U(4) = -8.9*SINN + 0.7*SIN2N
         U(5) = -2.1*SINN
         U(6) = U(5)
         U(8) = -17.7*SINN + 0.7*SIN2N
         U(9) = ATAN2(0.189*SINN - 0.0058*SIN2N,
     *                1.0 + 0.189*COSN - 0.0058*SIN2N)/RAD
         U(10) = U(9)
         U(11) = U(9)
         U(12) = ATAN2( 0.185*SINN, 1.0 + 0.185*COSN)/RAD
         U(13) = ATAN2(-0.201*SINN, 1.0 + 0.201*COSN)/RAD
         U(14) = ATAN2(-0.221*SINN, 1.0 + 0.221*COSN)/RAD
         U(18) = ATAN2(-0.198*SINN, 1.0 + 0.198*COSN)/RAD
         U(19) = ATAN2(-0.640*SINN - 0.134*SIN2N,
     *                 1.0 + 0.640*COSN + 0.134*COS2N)/RAD
         U(20) = ATAN2(-0.0373*SINN, 1.0 - 0.0373*COSN)/RAD
         U(21) = U(20)
         U(22) = U(20)
         U(24) = U(20)
         U(25) = ATAN2(-0.441*SINN, 1.0 + 0.441*COSN)/RAD
      ENDIF
*
*     sum over all tides
*     ------------------
      SUM = 0.D0
      DO 100 I=1,NT
         H1 = H12(1,I)
         H2 = H12(2,I)
         CHIU = (ARG(I) + U(I))*RAD
         SUM = SUM + H1*F(I)*COS(CHIU) + H2*F(I)*SIN(CHIU)
 100  CONTINUE
*
      TIDE = SUM
      RETURN
      END
*
*----------------------------------------------------------------------
      SUBROUTINE TLOAD( ARRAY,MX,MY, SELECT, LU, ZNULL,
     *                  NX,NY,LATMIN,LATMAX,LONMIN,LONMAX,
     *                  NGRIDS )
*
*  Loads ARRAY with the data from the input files of harmonic constants.
*
*  ARRAY  - holds inphase & quadrature components for all tides.
*  ZNULL  - desired value of NULL data.
*  LU     - tide coeffs will be read on fortran unit LU.
*  SELECT - 8-element logical array. NOT USED WITH PERTH2.
*
*  All remaining arguments describe the tidal grids, except:
*  NGRIDS - returns the total number of constituents loaded.
*
*
      REAL       ARRAY(MX,MY,2,*)
      REAL       LATMIN,LATMAX,LONMIN,LONMAX
      PARAMETER (NT=8)
      LOGICAL    SELECT(NT)
      PARAMETER (D2R=1.745329 E-2,  MXX=720,MXY=361)
      DIMENSION  AMP(MXX,MXY), PHA(MXX,MXY)
      CHARACTER  TITLE(2,8)*160
*
*     Read grids of ascii data
*     ------------------------
      NGRIDS = 0
      DO 300 IT=1,NT
         IF (SELECT(IT)) THEN
            CALL UTCSRI( AMP,MXX,MXY,NX,NY,TITLE(1,IT),LU,
     *                   LATMIN,LATMAX,LONMIN,LONMAX,UNDEF )
            CALL UTCSRI( PHA,MXX,MXY,NX,NY,TITLE(2,IT),LU,
     *                   LATMIN,LATMAX,LONMIN,LONMAX,UNDEF )
            NGRIDS = NGRIDS + 1
            IF (NX.GT.MX .OR. NY.GT.MY) THEN
               WRITE(6,5) MX,MY
    5          FORMAT('Size of ARRAY(',I4,',',I4,') must be increased.')
               STOP
            ENDIF
            DO 200 J=1,NY
            DO 200 I=1,NX
               IF (PHA(I,J).EQ.UNDEF .OR. AMP(I,J).EQ.UNDEF) THEN
                  ARRAY(I,J,1,NGRIDS) = ZNULL
                  ARRAY(I,J,2,NGRIDS) = ZNULL
                  ARRAY(I,J,1,1     ) = ZNULL
               ELSE
                  ARRAY(I,J,1,NGRIDS) = AMP(I,J) * COS(PHA(I,J)*D2R)
                  ARRAY(I,J,2,NGRIDS) = AMP(I,J) * SIN(PHA(I,J)*D2R)
               ENDIF
  200       CONTINUE
         ENDIF
  300 CONTINUE
  400 CONTINUE
*
*     Check for possible user setup error
*     -----------------------------------
      CALL CHECK1( SELECT, TITLE )
      RETURN
      END
*
      SUBROUTINE CHECK1( SELECT, TITLE )
*
*  Under the assumption that the TITLEs of the gridded input data
*  have the names of the tidal constituents somewhere in them,
*  this routine makes an elementary check for a user setup error.
*  If the tidal constituent name is not detected, it is assumed
*  that an error has occurred and the program is stopped.
*
      PARAMETER     (NT=8)
      LOGICAL        SELECT(NT)
      CHARACTER*160  TITLE(2,NT)
      INTRINSIC      INDEX
      LOGICAL        STOP
      CHARACTER*2    NAMES(NT)
      DATA   NAMES/'Q1','O1','P1','K1','N2','M2','S2','K2'/
*
      STOP = .FALSE.
      DO 10 I=1,NT
         IF (.NOT.SELECT(I)) GO TO 10
         IF (INDEX(TITLE(1,I),NAMES(I)).EQ.0) STOP = .TRUE.
         IF (INDEX(TITLE(2,I),NAMES(I)).EQ.0) STOP = .TRUE.
 10   CONTINUE
*
      IF (STOP) THEN
         WRITE(6,20)
 20      FORMAT(/' Subroutine CHECK1 has detected a likely error in'
     *        ' the input data.'/' Possible causes:'/
     *       4X,'The grids were read in the wrong order.'/
     *       4X,'A grid was input for a tide not SELECTed.'/
     *      ' Expected:',5X,'But read:')
         DO 30 I=1,NT
            IF (.NOT.SELECT(I)) GO TO 30
            WRITE(6,35) NAMES(I),TITLE(1,I)(1:60),
     *                           TITLE(2,I)(1:60)
 30      CONTINUE
 35      FORMAT(5X,A2,5X,A60/12X,A60)
         STOP 1
      ENDIF
      RETURN
      END
*===================================================================
      SUBROUTINE ASTROL( TIME, SHPN )
*
*  Computes the basic astronomical mean longitudes  s, h, p, N.
*  Note N is not N', i.e. N is decreasing with time.
*  These formulae are for the period 1990 - 2010, and were derived
*  by David Cartwright (personal comm., Nov. 1990).
*  TIME is UTC in decimal MJD.
*  All longitudes returned in degrees.
*  R. D. Ray    Dec. 1990
*
*  Non-vectorized version.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  SHPN(4)
      PARAMETER  (CIRCLE=360.0D0)
*
      T = TIME - 51544.4993D0
*
*     mean longitude of moon
*     ----------------------
      SHPN(1) = 218.3164D0 + 13.17639648D0 * T
*
*     mean longitude of sun
*     ---------------------
      SHPN(2) = 280.4661D0 +  0.98564736D0 * T
*
*     mean longitude of lunar perigee
*     -------------------------------
      SHPN(3) =  83.3535D0 +  0.11140353D0 * T
*
*     mean longitude of ascending lunar node
*     --------------------------------------
      SHPN(4) = 125.0445D0 -  0.05295377D0 * T

      SHPN(1) = MOD(SHPN(1),CIRCLE)
      SHPN(2) = MOD(SHPN(2),CIRCLE)
      SHPN(3) = MOD(SHPN(3),CIRCLE)
      SHPN(4) = MOD(SHPN(4),CIRCLE)

      IF (SHPN(1).LT.0.D0) SHPN(1) = SHPN(1) + CIRCLE
      IF (SHPN(2).LT.0.D0) SHPN(2) = SHPN(2) + CIRCLE
      IF (SHPN(3).LT.0.D0) SHPN(3) = SHPN(3) + CIRCLE
      IF (SHPN(4).LT.0.D0) SHPN(4) = SHPN(4) + CIRCLE
      RETURN
      END
*
*--------------------------------------------------------------------
      SUBROUTINE GRSINT( GRID,ND1,ND2,NGRIDS,NX,NY,UNDEF,
     *                   LATMIN,LATMAX,LONMIN,LONMAX,WRAP,
     *                   DLAT,DLON,VAL,ISDATA )
*
*  Interpolates a value from a grid of data at the desired location.
*  Interpolation is bilinear.
*  First 12 arguments above all describe the grid.
*  WRAP is T if grid allows wrap-around in longitude.
*  DLAT,DLON - location of desired position for interpolation.
*  VAL - returns interpolated value.
*  ISDATA is returned F if no valid data at position (DLAT,DLON).
*
*  R. Ray   3/13/91
*
*  Revised 6/3/93 to allow multiple grids to be handled.
*     VAL is now returned as an array; GRID is changed to 3 dimensions.
*     ISDATA is still a scalar (all grids assumed to have same nulls).
*
      REAL       LATMIN,LATMAX,LONMIN,LONMAX
      DIMENSION  GRID(ND1,ND2,*), VAL(*)
      LOGICAL    ISDATA,WRAP
      ISDATA = .TRUE.
      DX = (LONMAX - LONMIN)/REAL(NX - 1)
      DY = (LATMAX - LATMIN)/REAL(NY - 1)

*     compute indices for desired position
*     ------------------------------------
      JLAT1 = INT((DLAT - LATMIN)/DY) + 1
      JLAT2 = JLAT1 + 1
      IF (JLAT1.LT.1 .OR. JLAT2.GT.NY) THEN
         ISDATA = .FALSE.
         RETURN
      ENDIF
      XLON = DLON
      IF (WRAP.AND.XLON.LT.LONMIN) XLON = XLON + 360.0
      IF (WRAP.AND.XLON.GT.LONMAX) XLON = XLON - 360.0
      ILON1 = INT((XLON - LONMIN)/DX) + 1
      ILON2 = ILON1 + 1
      IF (ILON2.GT.NX) THEN
         IF (WRAP) THEN
            ILON2 = 1
         ELSE
            ISDATA = .FALSE.
            RETURN
         END IF
      ENDIF
      IF (ILON1.LT.1 .OR. ILON1.GT.NX) STOP 301 ! should never happen.

      IF (GRID(ILON1,JLAT1,1).EQ.UNDEF .AND.
     *    GRID(ILON2,JLAT1,1).EQ.UNDEF .AND.
     *    GRID(ILON1,JLAT2,1).EQ.UNDEF .AND.
     *    GRID(ILON2,JLAT2,1).EQ.UNDEF)      THEN
         ISDATA = .FALSE.
         RETURN
      ENDIF
      W1 = 0.0
      W2 = 0.0
      WX1 = (DX - (XLON - REAL(ILON1-1)*DX - LONMIN))/DX
      WX2 = 1.0 - WX1
      WY1 = (DY - (DLAT - REAL(JLAT1-1)*DY - LATMIN))/DY
      WY2 = 1.0 - WY1
*  Interpolation weights:
*  W1,W2,W3,W4 are for northwest,northeast,southeast,southwest corners.
      W1 = WX1*WY2
      W2 = WX2*WY2
      W3 = WX2*WY1
      W4 = WX1*WY1
      W = 0.0
      DO 10 I=1,NGRIDS
 10   VAL(I) = 0.0

*     get weights & interpolate
*     -------------------------
      IF (GRID(ILON1,JLAT1,1).NE.UNDEF) THEN
         W = W4
         DO 20 I=1,NGRIDS
 20      VAL(I) = W4*GRID(ILON1,JLAT1,I)
      ENDIF
      IF (GRID(ILON1,JLAT2,1).NE.UNDEF) THEN
         W = W + W1
         DO 30 I=1,NGRIDS
 30      VAL(I) = VAL(I) + W1*GRID(ILON1,JLAT2,I)
      ENDIF
      IF (GRID(ILON2,JLAT2,1).NE.UNDEF) THEN
         W = W + W2
         DO 40 I=1,NGRIDS
 40      VAL(I) = VAL(I) + W2*GRID(ILON2,JLAT2,I)
      ENDIF
      IF (GRID(ILON2,JLAT1,1).NE.UNDEF) THEN
         W = W + W3
         DO 50 I=1,NGRIDS
 50      VAL(I) = VAL(I) + W3*GRID(ILON2,JLAT1,I)
      ENDIF
      IF (W.GT.0.5) THEN
         DO 60 I=1,NGRIDS
 60      VAL(I) = VAL(I)/W
      ELSE
         ISDATA = .FALSE.
      ENDIF
      RETURN
      END
*===================================================================
      SUBROUTINE UTCSRI( G,NDIM1,NDIM2,NX,NY,TITLE,LU,
     *                   LATMIN,LATMAX,LONMIN,LONMAX,UNDEF )
*
*  Reads an ascii gridfile.
*  This format uses 80-byte card-image records in ascii.
*
*  Calling arguments:
*
*   G     - O - two-dimensional array to contain the gridded data.
*               The grid will be loaded in the following way:
*               G(i,j), with i going west to east and j south to north.
*               G(1,1) is thus the southwest corner of grid.
*               G(NX,NY) is the northeast corner.
*               User must make certain G is dimensioned to at least
*               the expected values of NX,NY.
*
*   NDIM1,2 - I - dimensions of G in calling program.
*               NDIM must not be less than the expected value of NX.
*
*   NX,NY - O - size of the grid G.
*
*   TITLE - O - 160-byte character-string title that describes data.
*
*   LU    - I - fortran unit to use for reading input data.
*
*   LATMIN,LATMAX,LONMIN,LONMAX - O - area limits for grid,
*               in decimal degrees.
*               Note: these are REAL variables.
*             The grid intervals are therefore:
*             DX = (LONMAX-LONMIN)/(NX-1)  and
*             DY = (LATMAX-LATMIN)/(NY-1).
*
*   UNDEF - O - value to denote a null or missing value in G.
*
*
*  Written by R. Ray      Aug. 1990
*  Modified 6/3/93 - added NDIM2 in calling arguments to test size of NY
*

      DIMENSION   G(NDIM1,NDIM2)
      REAL        LATMIN,LATMAX,LONMIN,LONMAX
      CHARACTER   TITLE*160, FORMAT*80

      READ(LU,1) TITLE(1:80)
      READ(LU,1) TITLE(81:160)
    1 FORMAT(A80)

      READ(LU,2) NY,NX
    2 FORMAT(16X,I5,16X,I5)
      IF (NX.GT.NDIM1 .OR. NY.GT.NDIM2) THEN
         WRITE(6,*) 'Increase dimensions in call to UTCSRI'
         STOP
      ENDIF
      READ(LU,3) LATMIN,LATMAX
      READ(LU,3) LONMIN,LONMAX
    3 FORMAT(16X,F9.0,16X,F9.0)
      READ(LU,3) UNDEF           ! 2 masks may be here
      READ(LU,1) FORMAT
      DO 4 J=1,NY
         READ(LU,FORMAT) (G(I,J),I=1,NX)
    4 CONTINUE
      RETURN
      END
