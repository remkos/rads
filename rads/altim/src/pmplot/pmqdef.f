**PMQDEF -- inquire projection type and scale
*+
      SUBROUTINE PMQDEF (TYPE, SCALE)
      CHARACTER*(*) TYPE
      REAL SCALE
*
* Query the type of geographical projection defined in PMDEF and the
* scale of the map defined in PMDEF or computed by PMWINDOW.
*
* Arguments:
*  TYPE   (output) : Character string quoting the projection type used.
*  SCALE  (output) : Scale (1:SCALE) of the map
*--
*  9-Jan-1991 - created [Remko Scharroo]
* 16-Jan-1991 - support conic and azimuthal projections.
* 30-Oct-1991 - support more projections.
* 19-Mar-1992 - Standardize PMPLOT.
*-----------------------------------------------------------------------
      INCLUDE 'pmplot.inc'
      TYPE=PTYPEC
      SCALE=PSCALE
      END
