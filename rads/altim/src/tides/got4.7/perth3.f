      SUBROUTINE PERTH3( DLAT,DLON,TIME,TIDE,ISDATA )
*
*  Name - PERTH                  Name derivation - PREdict Tidal Heights
*
*  Function -  to compute the ocean tidal height at a given time
*              and location from grids of harmonic constants.
*              Current version uses the 8 largest constituents in the
*              semidiurnal & diurnal bands, with other tides inferred,
*              plus the one radiational tide S1, plus optionally the
*              one compound tide M4.
*
*              (Long period tides are NOT computed by this routine.)
*
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
*     All input tide files should be concatenated, in order of frequency,
*     e.g., do:
*         cat q1.d o1.d p1.d s1.d k1.d n2.d m2.d s2.d k2.d > fort.30  OR
*         cat q1.d o1.d p1.d s1.d k1.d n2.d m2.d s2.d k2.d m4.d > fort.30
*     (or the equivalent for your operating system).
*     Notice that M4 is added OPTIONALLY at the end of the concatenation.
*
*  Processing logic -
*     Tidal constants at desired location are found by bilinear
*     interpolation from input grid files.  Nodal corrections are applied to 
*     all lunar tides.  Sixteen minor tides are inferred from major tides.
*
*  File references -
*     Input tidal grids are read on first call via unit LU (set in DATA).
*
*  Important local variables -
*     ARRAY - holds inphase & quadrature tidal constants as follows:
*        ARRAY(i,j,k,l) where
*        i,j are longitude/latitude indices.
*        k = 1,2  for Hcos(G) or Hsin(G).
*        l = 1,...,10 (max) - for each tidal constituent.
*     SHPN - holds the mean astronomical longitudes of interest;
*        this array is equivalenced to S, H, P, & OMEGA
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
*     2.1    1/13/98     R Ray    Fixed bad argument to "pi1" tide.
*     2.2    5/26/99     R Ray    Minor change: use different name IO routine.
*     2.3    5/27/99     R Ray    Added facility to pass dsn to TLOAD.
*     3.0    6/29/07     R Ray    Added handling of S1 and optional M4.
*                                 Array sizes increased for future use.
*                                 Astronomical longitudes now computed by ASTRO5.
*                                 Nodal factors now updated at every new day.
*     3.1    9/12/08     R Ray    Fixed bug on SHPN array size; let PP vary.
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  NINETY
      PARAMETER (MX=1440, MY=721, NG=10, NT=28)
      LOGICAL    INIT,ISDATA,SELECT(NG),WRAP,ISDATA0
      REAL       ARRAY(MX,MY,2,NG)
      REAL       LATMIN,LATMAX,LONMIN,LONMAX,UNDEF,H12(2,NT)
      REAL       SINN,SIN2N,COSN,COS2N
      DIMENSION  SHPN(5), ARG(NT), F(NT), U(NT)
      EQUIVALENCE (SHPN(1),S),(SHPN(2),H),(SHPN(3),P),(SHPN(4),OMEGA)
      EQUIVALENCE (SHPN(5),PP)
      SAVE
      PARAMETER (UNDEF=99999.0, PI=3.141592654D0, RAD=PI/180.D0)
      PARAMETER (TWO=2.D0, THREE=3.D0, FOUR=4.D0)
      PARAMETER (FIFTEN=15.D0, THIRTY=30.D0, NINETY=90.D0)
      DATA       DLAT0,DLON0,TIME0,TIME1/4*999.d0/
      DATA       INIT/.TRUE./, LU/30/, SELECT/10*.TRUE./, NX/0/

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
*         move "new" tides to end of H12 array (thus avoids some recoding)
*         ----------------------------------------------------------------
          h12(1,27) = h12(1,4)
          h12(2,27) = h12(2,4)
          h12(1,28) = h12(1,10)
          h12(2,28) = h12(2,10)
          h12(1,4) = h12(1,5)
          h12(2,4) = h12(2,5)
          h12(1,5) = h12(1,6)
          h12(2,5) = h12(2,6)
          h12(1,6) = h12(1,7)
          h12(2,6) = h12(2,7)
          h12(1,7) = h12(1,8)
          h12(2,7) = h12(2,8)
          h12(1,8) = h12(1,9)
          h12(2,8) = h12(2,9)
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
         CALL ASTRO5( TIME, SHPN )

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
         ARG(27) = T1 + 180.D0                         ! S1 (Doodson's phase)
         ARG(28) = TWO*ARG(6)                          ! M4
      ENDIF
*
*     determine nodal corrections f and u    
*        Note: Update this code next iteration of model  -RDR
*     -----------------------------------
      IF (ABS(TIME-TIME1).GT.1.D0) THEN
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
         F(28) = F(6)*F(6)

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
         U(28) = TWO*U(6)
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
      block data for_tload
      character   dsn*60
      common /tid_dsn_blk/ dsn
      data dsn/' '/
      end
*
*
      SUBROUTINE TLOAD( ARRAY,MX,MY, SELECT, LU, ZNULL,
     *                  NX,NY,LATMIN,LATMAX,LONMIN,LONMAX,
     *                  NGRIDS )
*
*  Loads ARRAY with the data from the input files of harmonic constants.
*
*  ARRAY  - holds inphase & quadrature components for all tides.
*  ZNULL  - desired value of NULL data.
*  LU     - tide coeffs will be read on fortran unit LU.
*  SELECT - 8-element logical array.  SET ALL TRUE IN SUBROUTINE PERTH2.
*
*  All remaining arguments describe the tidal grids, except:
*  NGRIDS - returns the total number of constituents loaded.
*
*  NOTE: SELECT is no longer used in the driver routine.
*
*  Revised 5/27/99
*    - If dsn is defined (by calling programs), then input file is
*      opened to that file.  Otherwise, this routine assumes file is
*      already opened (or opened to fort.LU).
*  Revised 6/29/07 - increased dimensions MX,MY,NT
*
      REAL       ARRAY(MX,MY,2,*)
      REAL       LATMIN,LATMAX,LONMIN,LONMAX
      PARAMETER (NT=10)
      LOGICAL    SELECT(NT)
      PARAMETER (D2R=1.745329 E-2,  MXX=1440,MXY=721)
      DIMENSION  AMP(MXX,MXY), PHA(MXX,MXY)
      CHARACTER  TITLE(2,NT)*160
      character   dsn*60
      common /tid_dsn_blk/ dsn
*
      if (dsn(1:1).ne.' ') open( lu, file=dsn, status='old' )
*
*     Read grids of ascii data
*     ------------------------
      NGRIDS = 0
      DO 300 IT=1,NT
         IF (SELECT(IT)) THEN
            CALL GRDINP( AMP,MXX,MXY,NX,NY,TITLE(1,IT),LU,
     *                   LATMIN,LATMAX,LONMIN,LONMAX,UNDEF )
            CALL GRDINP( PHA,MXX,MXY,NX,NY,TITLE(2,IT),LU,
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
      if (dsn(1:1).ne.' ') close(lu)
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
*  Revised 6/16/00 - allow for lower-case names.
*  Revised 6/29/07 - added S1.
*
      PARAMETER     (NT=9)
      LOGICAL        SELECT(*)
      CHARACTER*160  TITLE(2,*)
      INTRINSIC      INDEX
      LOGICAL        STOP
      CHARACTER*2    NAMES(NT), NAMEL(NT)
      DATA   NAMES/'Q1','O1','P1','S1','K1','N2','M2','S2','K2'/
      DATA   NAMEL/'q1','o1','p1','s1','k1','n2','m2','s2','k2'/
*
      STOP = .FALSE.
      DO 10 I=1,NT
         IF (.NOT.SELECT(I)) GO TO 10
         IF (INDEX(TITLE(1,I),NAMES(I)).EQ.0 .AND.
     *       INDEX(TITLE(1,I),NAMEL(I)).EQ.0) STOP = .TRUE.
         IF (INDEX(TITLE(2,I),NAMES(I)).EQ.0 .AND.
     *       INDEX(TITLE(2,I),NAMEL(I)).EQ.0) STOP = .TRUE.
 10   CONTINUE
*
      IF (STOP) THEN
         WRITE(6,20)
 20      FORMAT(/' Subroutine CHECK1 has detected a likely error in',
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


      SUBROUTINE GRDINP( G,NDIM1,NDIM2,NX,NY,TITLE,LU,
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
*  Revised 6/29/07 - if EOF hit, then set grids to zero.
*

      DIMENSION   G(NDIM1,NDIM2)
      REAL        LATMIN,LATMAX,LONMIN,LONMAX
      CHARACTER   TITLE*160, FORMAT*80
      LOGICAL     EOF
      DATA EOF/.FALSE./

      IF (EOF) GO TO 101

      READ(LU,1,END=101) TITLE(1:80)
      READ(LU,1) TITLE(81:160)
    1 FORMAT(A80)

      READ(LU,2) NY,NX
    2 FORMAT(16X,I5,16X,I5)
      IF (NX.GT.NDIM1 .OR. NY.GT.NDIM2) THEN
         WRITE(6,*) 'Increase dimensions in call to GRDINP'
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

*  Code added for optionally missing grids.....
 101  EOF = .true.
      if (NX.EQ.0) THEN
         WRITE(6,*) 'End-of-file hit on 1st call to PERTH READ routine.'
         WRITE(6,*) 'Abnormal STOP.'
         STOP 1
      ENDIF
      TITLE = ' '
      DO 120 J=1,NY
      DO 120 I=1,NX
         G(I,J) = 0.0
 120  CONTINUE
      RETURN
      END




      SUBROUTINE ASTRO5( TIME, SHPNP )
*
*---------------------------------------------------------------------
*  Computes the 5 basic astronomical mean longitudes  s, h, p, N, p'.
*
*  Note N is not N', i.e. N is decreasing with time.
*
*  TIME is UTC in decimal Modified Julian Day (MJD).
*  All longitudes returned in degrees.
*
*  R. D. Ray, NASA/GSFC   August 2003
*
*  Most of the formulae for mean longitudes are extracted from 
*  Jean Meeus, Astronomical Algorithms, 2nd ed., 1998.  
*  Page numbers below refer to this book.
*
*  Note: This routine uses TIME in UT and does not distinguish between
*    the subtle differences of UTC, UT1, etc.  This is more than adequate
*    for the calculation of these arguments, especially in tidal studies.
*---------------------------------------------------------------------
*
      IMPLICIT NONE
      DOUBLE PRECISION TIME, SHPNP(5)

      DOUBLE PRECISION TJD,T,CIRCLE,DEL,TJLAST,D
      PARAMETER       (CIRCLE=360.0D0)
      DOUBLE PRECISION DELTAT
      EXTERNAL         DELTAT
      SAVE             DEL,TJLAST
      DATA TJLAST/-1/

*     Convert to Julian Day and to Ephemeris Time
*     -------------------------------------------
      TJD = TIME + 2400000.5D0
      IF (ABS(TJD-TJLAST).GT.100.D0) THEN
         DEL = DELTAT( TJD )/86400.D0
         TJLAST = TJD
      ENDIF
      TJD = TJD + DEL

*     Compute time argument in centuries relative to J2000
*     ----------------------------------------------------
      T = ( TJD - 2451545.d0 )/36525.d0
*
*
*     mean longitude of moon (p.338)
*     ------------------------------
      SHPNP(1) = (((-1.53388d-8*T + 1.855835d-6)*T - 1.5786d-3)*T +
     *            481267.88123421d0)*T + 218.3164477d0
*
*     mean elongation of moon (p.338)
*     -------------------------------
      D = (((-8.8445d-9*T + 1.83195d-6)*T - 1.8819d-3)*T +
     *          445267.1114034d0)*T + 297.8501921d0
*
*     mean longitude of sun
*     ---------------------
      SHPNP(2) = SHPNP(1) - D
*
*     mean longitude of lunar perigee (p.343)
*     ---------------------------------------
      SHPNP(3) = ((-1.249172d-5*T - 1.032d-2)*T + 4069.0137287d0)*T +
     *            83.3532465d0
*
*     mean longitude of ascending lunar node (p.144)
*     ----------------------------------------------
      SHPNP(4) = ((2.22222d-6*T + 2.0708d-3)*T - 1934.136261d0)*T +
     *            125.04452d0
*
*     mean longitude of solar perigee (Simon et al., 1994)
*     ----------------------------------------------------
      SHPNP(5) = 282.94d0 + 1.7192d0 * T

      SHPNP(1) = MOD(SHPNP(1),CIRCLE)
      SHPNP(2) = MOD(SHPNP(2),CIRCLE)
      SHPNP(3) = MOD(SHPNP(3),CIRCLE)
      SHPNP(4) = MOD(SHPNP(4),CIRCLE)

      IF (SHPNP(1).LT.0.D0) SHPNP(1) = SHPNP(1) + CIRCLE
      IF (SHPNP(2).LT.0.D0) SHPNP(2) = SHPNP(2) + CIRCLE
      IF (SHPNP(3).LT.0.D0) SHPNP(3) = SHPNP(3) + CIRCLE
      IF (SHPNP(4).LT.0.D0) SHPNP(4) = SHPNP(4) + CIRCLE
      RETURN
      END







      double precision function deltat( tjd )

c..author: F X Timmes, University of Chicago

c..Slightly revised by R Ray, GSFC, Aug 2003.  Also updated table.

c..this routine computes the difference between 
c..universal time and dynamical time.

c..input:
c..tjd  = UT julian day number

c..output:
c  deltat = ET - UT in seconds


      implicit none
      save

c..declare the pass
      double precision tjd,tjde,secdif


c..for storing the table
      integer          tabstart,tabend,tabsize
      parameter       (tabstart = 1620,
     1                 tabend   = 2004,
     2                 tabsize  = tabend - tabstart + 1)
      double precision dt(tabsize),years(tabsize)


c..local variables
      integer          i,j,iy,im,id,iflag,mp,iat
      parameter        (mp = 2)
      double precision hr,b,xx,dy

c..parameter mp above determines the order of the interpolatant
c..when interpolating in the table. mp=2=linear mp=3=quadratic and so on.
c..don't be too greedy about the order of the interpolant since
c..for many years the data is flat, and trying to fit anything other
c..that a line will produce unwanted oscillations. thus i choose mp=2.
   


c..these tables of observed and extrapolated data are
c..are from the us naval observatory ftp://maia.usno.navy.mil/ser7/
c..years 1620 to 1710
      data  (dt(i), i=1,90) /
     1    124.00d0, 119.00d0, 115.00d0, 110.00d0, 106.00d0, 102.00d0, 
     2     98.00d0,  95.00d0,  91.00d0,  88.00d0,  85.00d0,  82.00d0,  
     3     79.00d0,  77.00d0,  74.00d0,  72.00d0,  70.00d0,  67.00d0,  
     4     65.00d0,  63.00d0,  62.00d0,  60.00d0,  58.00d0,  57.00d0, 
     5     55.00d0,  54.00d0,  53.00d0,  51.00d0,  50.00d0,  49.00d0, 
     6     48.00d0,  47.00d0,  46.00d0,  45.00d0,  44.00d0,  43.00d0,  
     7     42.00d0,  41.00d0,  40.00d0,  38.00d0,  37.00d0,  36.00d0, 
     8     35.00d0,  34.00d0,  33.00d0,  32.00d0,  31.00d0,  30.00d0, 
     9     28.00d0,  27.00d0,  26.00d0,  25.00d0,  24.00d0,  23.00d0, 
     &     22.00d0,  21.00d0,  20.00d0,  19.00d0,  18.00d0,  17.00d0,
     1     16.00d0,  15.00d0,  14.00d0,  14.00d0,  13.00d0,  12.00d0, 
     2     12.00d0,  11.00d0,  11.00d0,  10.00d0,  10.00d0,  10.00d0, 
     3      9.00d0,   9.00d0,   9.00d0,   9.00d0,   9.00d0,   9.00d0, 
     4      9.00d0,   9.00d0,   9.00d0,   9.00d0,   9.00d0,   9.00d0, 
     5      9.00d0,   9.00d0,   9.00d0,   9.00d0,  10.00d0,  10.00d0/


c..years 1711 to 1799
      data  (dt(i), i=91, 180) /
     1     10.00d0,  10.00d0,  10.00d0,  10.00d0,  10.00d0,  10.00d0, 
     2     10.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0, 
     3     11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0, 
     4     11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0, 
     5     12.00d0,  12.00d0,  12.00d0,  12.00d0,  12.00d0,  12.00d0,
     6     12.00d0,  12.00d0,  12.00d0,  12.00d0,  13.00d0,  13.00d0, 
     7     13.00d0,  13.00d0,  13.00d0,  13.00d0,  13.00d0,  14.00d0,  
     8     14.00d0,  14.00d0,  14.00d0,  14.00d0,  14.00d0,  14.00d0,  
     9     15.00d0,  15.00d0,  15.00d0,  15.00d0,  15.00d0,  15.00d0,  
     &     15.00d0,  16.00d0,  16.00d0,  16.00d0,  16.00d0,  16.00d0,
     1     16.00d0,  16.00d0,  16.00d0,  16.00d0,  16.00d0,  17.00d0,  
     2     17.00d0,  17.00d0,  17.00d0,  17.00d0,  17.00d0,  17.00d0,  
     3     17.00d0,  17.00d0,  17.00d0,  17.00d0,  17.00d0,  17.00d0,  
     4     17.00d0,  17.00d0,  17.00d0,  17.00d0,  16.00d0,  16.00d0,  
     5     16.00d0,  16.00d0,  15.00d0,  15.00d0,  14.00d0,  14.00d0/


c..years 1800 to 1890
      data  (dt(i), i=181, 270) /
     1     13.70d0,  13.40d0,  13.10d0,  12.90d0,  12.70d0,  12.60d0,  
     2     12.50d0,  12.50d0,  12.50d0,  12.50d0,  12.50d0,  12.50d0,  
     3     12.50d0,  12.50d0,  12.50d0,  12.50d0,  12.50d0,  12.40d0,  
     4     12.30d0,  12.20d0,  12.00d0,  11.70d0,  11.40d0,  11.10d0,  
     5     10.60d0,  10.20d0,   9.60d0,   9.10d0,   8.60d0,   8.00d0, 
     6      7.50d0,   7.00d0,   6.60d0,   6.30d0,   6.00d0,   5.80d0,
     7      5.70d0,   5.60d0,   5.60d0,   5.60d0,   5.70d0,   5.80d0,
     8      5.90d0,   6.10d0,   6.20d0,   6.30d0,   6.50d0,   6.60d0,  
     9      6.80d0,   6.90d0,   7.10d0,   7.20d0,   7.30d0,   7.40d0,  
     &      7.50d0,   7.60d0,   7.70d0,   7.70d0,   7.80d0,   7.80d0, 
     1      7.88d0,   7.82d0,   7.54d0,   6.97d0,   6.40d0,   6.02d0,  
     2      5.41d0,   4.10d0,   2.92d0,   1.82d0,   1.61d0,   0.10d0,  
     3     -1.02d0,  -1.28d0,  -2.69d0,  -3.24d0,  -3.64d0,  -4.54d0,  
     4     -4.71d0,  -5.11d0,  -5.40d0,  -5.42d0,  -5.20d0,  -5.46d0,  
     5     -5.46d0,  -5.79d0,  -5.63d0,  -5.64d0,  -5.80d0,  -5.66d0/


c..years 1890 to 1979
      data  (dt(i), i=271, 360) /
     1     -5.87d0,  -6.01d0,  -6.19d0,  -6.64d0,  -6.44d0,  -6.47d0,  
     2     -6.09d0,  -5.76d0,  -4.66d0,  -3.74d0,  -2.72d0,  -1.54d0,  
     3      -0.2d0,   1.24d0,   2.64d0,   3.86d0,   5.37d0,   6.14d0,  
     4      7.75d0,   9.13d0,  10.46d0,  11.53d0,  13.36d0,  14.65d0,  
     5     16.01d0,  17.20d0,  18.24d0,  19.06d0,  20.25d0,  20.95d0, 
     6     21.16d0,  22.25d0,  22.41d0,  23.03d0,  23.49d0,  23.62d0,  
     7     23.86d0,  24.49d0,  24.34d0,  24.08d0,  24.02d0,  24.00d0,  
     8     23.87d0,  23.95d0,  23.86d0,  23.93d0,  23.73d0,  23.92d0,  
     9     23.96d0,  24.02d0,  24.33d0,  24.83d0,  25.30d0,  25.70d0,  
     &     26.24d0,  26.77d0,  27.28d0,  27.78d0,  28.25d0,  28.71d0, 
     1     29.15d0,  29.57d0,  29.97d0,  30.36d0,  30.72d0,  31.07d0,  
     2     31.35d0,  31.68d0,  32.18d0,  32.68d0,  33.15d0,  33.59d0,  
     3     34.00d0,  34.47d0,  35.03d0,  35.73d0,  36.54d0,  37.43d0,  
     4     38.29d0,  39.20d0,  40.18d0,  41.17d0,  42.23d0,  43.37d0,  
     5     44.49d0,  45.48d0,  46.46d0,  47.52d0,  48.53d0,  49.59d0/

c..years 1980 to 2004 
      data  (dt(i), i=361, 385) /
     1     50.54d0,  51.38d0,  52.17d0,  52.96d0,  53.79d0,  54.34d0,  
     2     54.87d0,  55.32d0,  55.82d0,  56.30d0,  56.86d0,  57.57d0,  
     3     58.31d0,  59.12d0,  59.98d0,  60.78d0,  61.63d0,  62.29d0, 
     4     62.97d0,  63.47d0,  63.83d0,  64.09d0,  64.30d0,  64.50d0,  
     5     65.00d0/
c  NOTE:   Values thru 2002 ok, 2003 & above are extrapolated. -RDR

c..first time flag
      data iflag/0/



c..if this is the first time, initialize 
c..put the julain date at the first of the year in the array
c..years so that we can interpolate/extrapolate directly on
c..the given ut julian date

      if (iflag .eq. 0) then
       iflag = 1
  
       do i=1,tabsize
        j = tabstart + (i-1)
        call juldat2(j,1,1,0.0d0,xx)
        years(i) = xx
       end do

      endif


c..convert the given ut julian date to a ut calander date
      call caldat2(tjd,iy,im,id,hr)


c..if we are outside the table on the low end
c..use the stephenson and morrison expression 948 to 1600,
c..and the borkowski formula for earlier years

      if (iy .lt. tabstart) then
       if (iy .ge. 948) then
        b      = 0.01d0 * float(iy - 2000)
        secdif = b * (b * 23.58d0 + 100.3d0) + 101.6d0
       else
        b      = 0.01d0 * float(iy -2000) + 3.75d0
        secdif = 35.0d0 * b * b + 40.0d0
       end if


c..if we are outside the table on the high end
c..use a linear extrapolation into the future

      else if (iy .gt. tabend) then
       b      = float(iy - tabend)
       secdif = dt(tabsize) + b*(dt(tabsize) - dt(tabsize-1))



c..otherwise we are in the table
c..get the table location and interpolate

      else
       iat = iy - tabstart + 1
       iat = max(1, min(iat - mp/2 + 1, tabsize - mp + 1))
       call polint(years(iat),dt(iat),mp,tjd,secdif,dy)


c..the astronomical almanac table is corrected by adding the expression
c..      -0.000091 (ndot + 26)(year-1955)^2  seconds
c..to entries prior to 1955 (page K8), where ndot is the secular tidal 
c..term in the mean motion of the moon. entries after 1955 are referred 
c..to atomic time standards and are not affected by errors in lunar 
c..or planetary theory.  a value of ndot = -25.8 arcsec per century squared 
c..is the value used in jpl's de403 ephemeris, the earlier de200 ephemeris
c..used the value -23.8946. note for years below the table (less than 1620)
c..the time difference is not adjusted for small improvements in the 
c..current estimate of ndot because the formulas were derived from 
c..studies of ancient eclipses and other historical information, whose 
c..interpretation depends only partly on ndot.
c..here we make the ndot correction.

       if (iy .lt. 1955) then
        b = float(iy - 1955)
        secdif = secdif - 0.000091d0 * (-25.8d0 + 26.0d0)*b*b
       end if
      end if



c..add the difference to the ut julian date to get the dynamical julian date
c     tjde = tjd + secdif/86400.0d0

      deltat = secdif

      return
      end




      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      save

c..given arrays xa and ya of length n and a value x, this routine returns a 
c..value y and an error estimate dy. if p(x) is the polynomial of degree n-1
c..such that ya = p(xa) ya then the returned value is y = p(x) 

c..input:
c..xa(1:n) = array of x values
c..ya(1:n) = array of y values
c..n       = order of interpolant, 2=linear, 3=quadratic ...
c..x       = x value where interpolation is desired

c..output:
c..y       = interpolated y value
c..dy      = error esimate


c..declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=10)
      double precision xa(n),ya(n),x,y,dy,c(nmax),d(nmax),dif,dift,
     1                 ho,hp,w,den

c..find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
      enddo

c..first guess for y
      y = ya(ns)

c..for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do m=1,n-1
       do i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
       enddo

c..after each column is completed, decide which correction c or d, to add
c..to the accumulating value of y, that is, which path to take in the table
c..by forking up or down. ns is updated as we go to keep track of where we
c..are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
      enddo
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c..The next two routines are alternatives to the juldat and caldat routines
c..in the novas package, since the novas routines do not take care of the
c..switch from julian to gregorian calanders on 15oct1582. Hence,
c..the novas routines are of no use with the long jpl ephemris de406.

c..Input calendar dates 1582-oct-15 and after are taken to be expressed 
c..in the gregorian calendar system. Prior dates are assumed to be in the 
c..julian calendar system.

c..Historically, not all regions switched calendars at the same time 
c..(or even in the same century). Thus, the user must be aware of
c..which calendar was in effect for a particular historical record. 
c..It should not be assumed this system's calendar automatically
c..correlates with a date from an arbitrary historical document. 

c..Here is the progression near the calendar switch point: 
c
c       calendar type    calendar date   julian day number
c       -------------    -------------   -----------------
c        julian           1582-oct-03        2299158.5
c        julian           1582-oct-04        2299159.5 --->
c         (skipped)      "1582-oct-05"                    |
c         (skipped)      "1582-oct-06"                    |
c         (skipped)      "1582-oct-07"                    |
c         (skipped)      "1582-oct-08"                    |
c         (skipped)      "1582-oct-09"                    |
c         (skipped)      "1582-oct-10"                    |
c         (skipped)      "1582-oct-11"                    |
c         (skipped)      "1582-oct-12"                    |
c         (skipped)      "1582-oct-13"                    |
c         (skipped)      "1582-oct-14"                    |
c        gregorian        1582-oct-15        2299160.5 <---
c        gregorian        1582-oct-16        2299161.5
c        gregorian        1582-oct-17        2299162.5


c..in this system there are zero and negative years. 
c..the progression is as follows: 
c
c    julian day number       labeling-convention
c      (jan 1 00:00)       BC/AD      arithmetical 
c    -----------------     -----      ------------
c    1720327.5              3BC           -2
c    1720692.5              2BC           -1
c    1721057.5              1BC            0
c    1721423.5              1AD            1
c    1721788.5              2AD            2


c  Author: F. X. Timmes, U. Chicago



      subroutine juldat2(iy,im,id,rh,tjd)
      implicit none
      save

c..converts a calander date to a julian date.
c..input time value can be in any ut-like time scale (utc, ut1, tt, etc.) 
c..and the output julian date will have same basis.

c..the astronomical calendar is used. thus the year before 1 AD is 0, 
c..the year before that is 1 BC. the change to the gregorian calander 
c..on oct 15, 1582 is accounted for.

c..input:
c..iy = integer year
c..im = integer month
c..id = integer day
c..rh = number of hours past midnight

c..output:
c..tjd = julian date


c..declare
      integer          im,iy,id,igreg,i,jd
      parameter        (igreg = 15+31*(10+12*1582))
      double precision tjd,xy,xm,xa,xb,rh

      if (im .gt. 2) then
       xy = iy
       xm = im + 1
      else
       xy = iy - 1
       xm = im + 13
      end if
      i = id + 31 * (im + 12 * iy)
      if (i .ge. igreg) then
       xa = int(0.01d0 * xy)
       xb = 2.0d0 - xa + int(0.25d0 * xa)
      else
       xb = 0.0d0
      end if
      jd  = int(365.25d0*(xy + 4716.0d0)) + int(30.6001d0*xm) + id
      tjd = dble(jd) + xb - 1524.5d0 + rh/24.d0
      return
      end




      subroutine caldat2(tjd,iy,im,id,rh)
      implicit none
      save

c..converts a julian date to a calander date.
c..input time value can be in any ut-like time scale (utc, ut1, tt, etc.) 
c..and the output calander date will have same basis.

c..the astronomical calendar is used. thus the year before 1 ad is 0, 
c..the year before that is 1 bc. the change to the gregorian calander 
c..on oct 15, 1582 is accounted for.

c..input:
c..tjd = julian date
c..iy = integer year
c..im = integer month
c..id = integer day
c..rh = number of hours past midnight

c..output:
c..iy = integer year
c..im = integer month
c..id = integer day
c..rh = number of hours past midnight


c..declare
      integer          id,im,iy,igreg
      parameter        (igreg = 2299161)
      double precision tjd,rh,x1,z,f,x2,xa,xb,xc,xd,xe,rd,c1,c2,c3
      parameter        (c1 = 1.0d0/36524.25d0,
     1                  c2 = 1.0d0/365.25d0, 
     2                  c3 = 1.0d0/30.6001d0) 

      x1 = tjd + 0.5d0
      z  = int(x1)
      f  = x1 - z

      if (x1 .ge. igreg) then
       x2 = int((x1-1867216.25d0)*c1)
       xa = z + 1.0d0 +x2 - int(0.25d0 * x2)
      else
       xa = z
      end if
      xb = xa + 1524.0d0
      xc = int((xb - 122.1d0)*c2)
      xd = int(365.25d0*xc)
      xe = int((xb-xd)*c3)

      rd = xb - xd - int(30.6001d0*xe) + f
      id = rd
      rh = 24.0d0*(rd - dble(id)) 
      im = xe - 1
      if (im .gt. 12) im = im - 12
      iy = xc - 4715
      if (im .gt. 2) iy = iy - 1
      return
      end




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
*  Revised 1/10/95 to test longitude limits when wrap=.false.
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
      if (.not.wrap .and. (ilon1.lt.1 .or. ilon2.gt.nx)) then
         isdata = .false.
         return
      endif
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
