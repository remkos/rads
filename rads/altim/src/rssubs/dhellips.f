**DHELLIPS -- Compute height difference between ellipsoids
*+
      FUNCTION DHELLIPS (CONV, LAT)
      REAL*8   DHELLIPS, LAT
      INTEGER  CONV

* Compute height difference between WGS84, GEOSAT and TOPEX ellipsoids.
* The input is LAT in degrees, the returned function value is DHELLIPS
* in metres (the result is always non-negative). CONV specifies the
* conversion (see below).
*
* The ellipsoids are defined by their equatorial radius and inverse
* flattening as shown in the table below. Also indicated are which
* products use which ellipsoids.
*
* Ellipsoid  Eq. radius (m)  Inv. flattening (-)  Products
* WGS84      6378137.0       298.257223563        ESA ERS and Envisat
* GEOSAT     6378137.0       298.257              GEOSAT T2 GDRs,
*                                                 DEOS ERS orbits
* TOPEX      6378136.3       298.257              AVISO T/P and Jason,
*                                                 GEOSAT J3 GDRs, RADS
*
* The conversion is linearised in sin(lat)**2, changes in latitude
* are ignored.
*
* Input Arguments:
*  CONV     : Conversion indicator
*             CONV=1 : Height of WGS84 over TOPEX ellipsoid
*             CONV=2 : Height of WGS84 over GEOSAT ellipsoid
*             CONV=3 : Height of GEOSAT over TOPEX ellipsoid
*             Other  : Zero height difference is returned
*  LAT      : Latitude (degrees)
*
* Function value (output):
*  DHELLIPS : Height difference between ellipsoids (m)
*
*
* Example 1: Convert DEOS orbital altitudes (provided relative to the
* GEOSAT ellipsoid) to the WGS84 reference ellipsoid.
*
*     CALL GETORB (TIME, DLAT, DLON, HORBIT, PATH, VERBOSE)
*     HORBIT = HORBIT - DHELLIPS (2, DLAT)
*
* Example 2: Convert heights in RA2 data products (WGS84 reference)
* to TOPEX ellipsoid reference.
*
*     IDH = NINT(DHELLIPS(1,LAT/1D6)*1D3)
*     ALT_COG_ELLIP = ALT_COG_ELLIP + IDH
*     GEOID_HT = GEOID_HT + IDH
*     M_SEA_SURF_HT = M_SEA_SURF_HT + IDH
*-
* 22-Feb-2002 - Created by Remko Scharroo
* 25-Feb-2002 - Comments updated
*-----------------------------------------------------------------------
      real*8 ae_top,f_top,ae_geo,f_geo,ae_wgs,f_wgs,pi,rad
      parameter (ae_top=6378136.3d0,f_top=1d0/298.257d0,
     |           ae_geo=6378137.0d0,f_geo=f_top,
     |           ae_wgs=ae_geo,     f_wgs=1d0/298.257223563d0)
      parameter (pi=3.14159265358979d0,rad=pi/180d0)

      if (conv.eq.1) then	! WGS84 - TOPEX
         dhellips=(ae_wgs-ae_top)+
     |            (ae_top*f_top-ae_wgs*f_wgs)*sin(lat*rad)**2
      else if (conv.eq.2) then	! WGS84 - GEOSAT
         dhellips=ae_geo*(f_geo-f_wgs)*sin(lat*rad)**2
      else if (conv.eq.3) then	! GEOSAT - TOPEX
         dhellips=(ae_geo-ae_top)*(1d0-f_geo*sin(lat*rad)**2)
      else
         dhellips=0d0
      endif
      end
