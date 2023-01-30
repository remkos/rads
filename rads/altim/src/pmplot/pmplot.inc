*-----------------------------------------------------------------------
* PMPLOT: common block definition for geographical projection routines
*-----------------------------------------------------------------------
* General parameters:
*
* PI      pi.
* QPI     quarter pi (pi/4).
* RAD     conversion degrees to rad (=pi/180).
* RAD2    = RAD/2
* RMEAN   mean earth radius in kilometers.
*
* Projection parameters:
* PVERSION PMPLOT Version number
* PMOPEN  identifier for the use of PMPLOT (2=PMDEF passed, 4=PMWINDOW
*	  passed)
* PTYPE   projection identifier from PMDEF.
* PTYPEC  projection type (character).
* PSCALE  scale of the map.
* PPARA1  projection parameter 1 (set in PMDEF)
* PPARA2  projection parameter 2 (set in PMDEF)
* FACTK   factor k (degrees or unity).
* FACTL   factor l (1/degrees or unity).
* RVIEW   height of observer above CM or maximum viewing range.
* AZMTAL  azimuthal projection selected.
* SINB0   = SIN(LATCEN*RAD)
* COSB0   = COS(LATCEN*RAD)
*
* Box parameters:
* XCURVE  specifies whether the parallels are curved.
* YCURVE  specifies whether the meridians are curved.
*
* Map boundaries:
* LONMIN  minimum longitude.
* LONCEN  center longitude.
* LONMAX  maximum longitude.
* LATMIN  minimum latitude.
* LATCEN  center latitude.
* LATMAX  maximum latitude.
* LONUMIN minimum longitude specified by user
* LONUMAX maximum longitude specified by user
* LATUMIN minimum latitude specified by user
* LATUMAX maximum latitude specified by user
*--
*  9-Jan-1991 - created [Remko Scharroo]
* 14-Jan-1991 - include conic and azimuthal projections.
* 19-Mar-1992 - add some parameters. Standardize PMPLOT.
*  6-Apr-1993 - add version number
*-----------------------------------------------------------------------
      INTEGER   PMOPEN,PTYPE
      CHARACTER PTYPEC*80
      LOGICAL   XCURVE,YCURVE,AZMTAL
      REAL      PI,QPI,RAD,RAD2,RMEAN,PVERSION
      REAL      PSCALE,PPARA1,PPARA2,FACTK,FACTL,RVIEW,SINB0,COSB0
      REAL      LONMIN,LONCEN,LONMAX,LATMIN,LATCEN,LATMAX
      REAL      LONUMIN,LONUMAX,LATUMIN,LATUMAX
*
      COMMON /PMPLT/ PI,QPI,RAD,RAD2,RMEAN,PVERSION,PMOPEN,
     .PTYPE,PTYPEC,PSCALE,PPARA1,PPARA2,FACTK,FACTL,RVIEW,SINB0,COSB0,
     .XCURVE,YCURVE,AZMTAL,
     .LONMIN,LONCEN,LONMAX,LATMIN,LATCEN,LATMAX,
     .LONUMIN,LONUMAX,LATUMIN,LATUMAX
      SAVE    /PMPLT/
