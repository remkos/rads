      SUBROUTINE LPEQMT( TS, DLAT, TLP )
*
*  Function - computes the long-period equilibrium ocean tides.
*
*  Arguments -
*     name      type  I/O               description
*     ----      ----  ---               -----------
*     TS         D     I   modified Julian day, in seconds, denoting
*                          time at which tide is to be computed.
*
*     DLAT       D     I   latitude in degrees (positive north) for the
*                          position at which tide is computed.
*
*     TLP        D     O   computed long-period tide, in centimeters.
*
*
*  Processing logic -
*     Fifteen tidal spectral lines from the Cartwright-Tayler-Edden
*     tables are summed over to compute the long-period tide.
*
*  Technical references -
*     Cartwright & Tayler, Geophys. J. R.A.S., 23, 45, 1971.
*     Cartwright & Edden, Geophys. J. R.A.S., 33, 253, 1973.
*
*  History -
*    version   date       programmer         change description
*    -------   ----       ----------         ------------------
*      1.0   11/27/90     D. Cartwright    Documented by R. Ray
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  PHC(4),DPD(4),SHPN(4)
      SAVE
*
      PARAMETER (RAD= 0.017453292519943 D0)
      PARAMETER (PSOL=283.D0*RAD)
      DATA PHC/290.21D0,280.12D0,274.35D0,343.51D0/,
     *     DPD/13.1763965D0,0.9856473D0,0.1114041D0,0.0529539D0/,
     *     TC/40431744.D2/
     
      INTEGER N
*
*     Compute 4 principal mean longitudes in radians at time TD
*     ---------------------------------------------------------
      TD = (TS - TC)/864.D2
      DO 1 N=1,4
         PH = PHC(N) + TD*DPD(N)
         SHPN(N) = RAD*DMOD(PH,36.D1)
    1 CONTINUE
*
*     Assemble long-period tide potential from 15 CTE terms > 1 mm.
*     Nodal term is included but not the constant term.
*     -------------------------------------------------
      ZLP = 2.79D0 * DCOS( SHPN(4) )
     *     -0.49D0 * DCOS( SHPN(2) - PSOL )
     *     -3.10D0 * DCOS( 2.D0*SHPN(2) )
      PH = SHPN(1)
      ZLP = ZLP - 0.67D0 * DCOS( PH - 2.D0*SHPN(2) + SHPN(3) )
     *          -(3.52D0 - 0.46D0*DCOS(SHPN(4))) * DCOS( PH - SHPN(3) )
      PH = PH + SHPN(1)
      ZLP = ZLP - 6.66D0 * DCOS( PH )
     *          - 2.76D0 * DCOS( PH + SHPN(4) )
     *          - 0.26D0 * DCOS( PH + 2.D0*SHPN(4) )
     *          - 0.58D0 * DCOS( PH - 2.D0*SHPN(2) )
     *          - 0.29D0 * DCOS( PH - 2.D0*SHPN(3) )
      PH = PH + SHPN(1)
      ZLP = ZLP - 1.27D0 * DCOS( PH - SHPN(3) )
     *          - 0.53D0 * DCOS( PH - SHPN(3) + SHPN(4) )
     *          - 0.24D0 * DCOS( PH - 2.D0*SHPN(2) + SHPN(3) )
*
*     Multiply by gamma_2 * sqrt(5/4 pi) * P20(lat)
*     ---------------------------------------------
      S = DSIN(DLAT*RAD)
      TLP = 0.437D0*ZLP*(1.5D0*S*S - 0.5D0)
      RETURN
      END
