**PMDEF -- set geographical projection type and scale
*+
      SUBROUTINE PMDEF (PROJECT, SCALE, PARA1, PARA2)
      INTEGER PROJECT
      REAL SCALE, PARA1, PARA2
*
* Define projection type and open possibility to use the PMPLOT routines.
* Also choose the map scale to be used and the standard or true-scale
* parallels.
* If you select any of the non-linear projection types, conversion
* between real-world coordinates (longitude,latitude) and map coordinates
* (x,y) is mandatory, and PMCONV must be used. For equi-rectangular projection
* it is optional.
* If you want to use any of the other PMPLOT routines, call PMDEF first.
*
* If you select a projection type less or equal to zero, standard PGPLOT
* routines will be called by PMWINDOW and PMBOX.
*
* Arguments:
*  PROJECT  (input) : Projection type:
*
*  [ <=0: Non-geographical projections:]
*   =  0, Adjusted window (PGWNAD) (same scale both direction)
*   = -1, Maximum window (PGSWIN) (unequal scales)
*
*  [1-10: Cylindrical projections:]
*   =  1, Equi-rectangular (linear).
*   =  2, Peters (non-linear).
*   =  3, Mercator (non-linear).
*   =  4, Miller (non-linear).
*   =  5, Gall's stereographic (non-linear).
*   =  6, ERS projection (non-linear).
*
* [11-20: Azimuthal projections (non-linear):]
*   = 11, Orthographic.
*   = 12, Perspective.
*   = 13, Azimuthal equal-area.
*   = 14, Azimuthal equi-distant.
*
* [21-30: Conic projections (non-linear):]
*   = 21, Ptolemy or conic equal-interval.
*   = 22, Kavraiskiy IV or conic equal-interval.
*
* [31-40: Miscellaneous projections (non-linear):]
*   = 31, Sinusoidal or Mercator equal-area.
*   = 32, Mollweide projection.
*   = 33, Bartholomew's `The Times' projection.
*   = 34, Tilted rectangular.
*
* [41-50: Polar projections (non-linear):]
*   = 41: Polar projection centred on the North Pole
*   = 42: Polar projection centred on the South Pole
*
* For all projections other than Azimuthal:
*  SCALE (input) : Scale of the map (1:SCALE) along the true-scale or standard
*                  parallel(s). If the plot does not have to be at any
*                  particular scale, specify SCALE=0.0.
*  PARA1 (input) : Latitude of the standard (or true-scale) parallel.
*  PARA2 (input) : (only for projections having two standard parallels; e.g.
*                  Kavraiskiy IV) Latitude of the second standard parallel.
*
* For Azimuthal projections:
*  SCALE (input) : Scale of the map (1:SCALE) at the center of the plot.
*                  If the plot does not have to be at any particular scale,
*                  specify SCALE=0.0.
*  PARA1, PARA2  : (not used).
*
* For tilted rectangular projection (34):
*  SCALE (input) : Scale of the map (1:SCALE) at the center of the plot.
*                  If the plot does not have to be at any particular scale,
*                  specify SCALE=0.0.
*  PARA1 (input) : Vertical scaling (unity).
*  PARA2 (input) : Tilt (degrees).
*--
*  9-Jan-1991 - created [Remko Scharroo]
* 14-Jan-1991 - include conic projections.
* 16-Jan-1991 - support azimuthal projections.
* 30-Oct-1991 - support more projections.
* 19-Mar-1992 - Standardize PMPLOT.
*  2-Apr-1993 - Tilted rectangular introduced. STPAR -> PARA, LATSP -> PPARA
* 30-Jun-1994 - Include polar projections.
* 11-Jul-1996 - Adjusted to PGPLOT 5.1
* 14-Aug-1997 - Add ERS projection
*-----------------------------------------------------------------------
      INCLUDE 'pgplot.inc'
      INCLUDE 'pmplot.inc'
      LOGICAL PGNOTO
      PVERSION=9607.0
*
      IF (PGNOTO('PMDEF')) RETURN
*
      QPI=ATAN(1.)
      PI=4*QPI
      RAD=PI/180
      RAD2=RAD/2
      RMEAN=6371.
*
      PMOPEN=2
      PTYPE=PROJECT
      PTYPEC='Use PMWINDOW first.'
      PSCALE=SCALE
      PPARA1=PARA1
      PPARA2=PARA2
      END
