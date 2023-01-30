**SETEARTH -- Specify parameters of the reference ellipsoid
*+
      SUBROUTINE SETEARTH (A_E, F_INV)
      REAL*8 A_E, F_INV

* Specify the equatorial radius (A_E) and inverse flattening (F_INV)
* of the reference ellipsoid to be used by GEOCEN, GEODET, GEOXYZ
* and XYZGEO.
*
* When this routine is not called the default values will be used:
* A_E = 6378137.0 meters
* F_INV = 298.257
*
* Input arguments:
*  A_E : Mean equatorial radius in meters
*  F_INV : Inverse flattening
*-
* 16-Jan-2002: Created by Remko Scharroo
*-----------------------------------------------------------------------
      include "earth.inc"
      AE=A_E
      FINV=F_INV
      end

**GETEARTH -- Inquire parameters of the reference ellipsoid
*+
      SUBROUTINE GETEARTH (A_E, F_INV)
      REAL*8 A_E, F_INV

* Get the equatorial radius (A_E) and inverse flattening (F_INV)
* of the reference ellipsoid to be used by GEOCEN, GEODET, GEOXYZ
* and XYZGEO.
*
* Output arguments:
*  A_E : Mean equatorial radius in meters
*  F_INV : Inverse flattening
*-
* 16-Jan-2002: Created by Remko Scharroo
*-----------------------------------------------------------------------
      include "earth.inc"
      A_E=AE
      F_INV=FINV
      end

      block data
      include "earth.inc"
      data ae, finv / 6378137d0, 298.257d0 /
      end
