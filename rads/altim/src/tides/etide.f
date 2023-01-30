**ETIDE -- Solid earth tidal elevation.
*+
      SUBROUTINE ETIDE (LAT, LON, TIME, ETSUN, ETMOON, ETFDEP, ETPERM)
      REAL*8 LAT, LON, TIME, ETSUN, ETMOON, ETFDEP, ETPERM
*
* Compute radial solid earth tidal displacement according to Love model.
* The output is the solid tidal (vertical) elevation on the ellipsoid,
* caused by the Sun and the Moon. Also the permanent tide and frequency
* dependent tide are given.
* According to the MERIT standards the Solid Earth tidal displacement
* does not include the permanent tide. So
*
*   ETIDE = ETSUN + ETMOON + ETFDEP - ETPERM
*
* The NOAA formulation of the solid earth tide does not include the
* frequency dependent tides. So
*
*   ETIDE = ETSUN + ETMOON - ETPERM
*
* GEODYN II, however, follows the IERS convention and has the permanent
* tide included. So
*
*   ETIDE = ETSUN + ETMOON + ETFDEP
*
* Arguments:
*  LAT     (input): Geodetic latitude (degrees) of the location at which the
*                   solid Earth tidal elevation must be computed.
*  LON     (input): East longitude (degrees) of this location.
*  TIME    (input): Time in MJD.
*  ETSUN  (output): Radial displacement (in meters) of the solid Earth at the
*                   given location due to gravitational attraction by the Sun.
*  ETMOON (output): Radial displacement (in meters) of the solid Earth at the
*                   given location due to gravitational attraction by the Moon.
*  ETFDEP (output): Frequency dependent solid Earth tide (in meters).
*  ETPERM (output): Time invariant radial solid earth tidal displacement (in
*                   meters) at the given latitude.
*-
* 22-Nov-1991: Created (Remko Scharroo).
*  9-Feb-1994: Revised.
*-----------------------------------------------------------------------
      include '../rssubs/etide.inc'
      include '../rssubs/math.inc'
      include '../rssubs/GRS80.inc'
      
      real*8 x(3),z(3),elem(6,2)
      real*8 a,anmtot,etide0,epharg,r,rlat,rlon

      elem(1,1)=amoon
      elem(2,1)=emoon
      elem(3,1)=imoon*rad
      elem(1,2)=asun
      elem(2,2)=esun
      elem(3,2)=isun*rad

      rlat=lat*rad
      rlon=(lon+epharg(1,time))*rad
      call geoxyz(rlat,rlon,0d0,x,r)
      
* Moon

      elem(4,1)=epharg(5,time)*rad
      a=epharg(4,time)*rad
      elem(5,1)=a-elem(4,1)
      elem(6,1)=anmtot(epharg(2,time)*rad-a,elem(2,1))
      call elevec(elem(1,1),z,-1d0)
      call rotate(1,-elem(3,2),z,z)
      etmoon=etide0(x,z,gmmoon,am,h2love)

* Sun

      elem(4,2)=0d0
      elem(5,2)=epharg(6,time)*rad
      elem(6,2)=anmtot(epharg(3,time)*rad-elem(5,2),elem(2,2))
      call elevec(elem(1,2),z,-1d0)
      etsun=etide0(x,z,gmsun,am,h2love)

* Permanent tide

      etperm=tperm*(1.5d0*sin(rlat)**2-0.5d0)

* Frequency dependent tide

      etfdep=tfdep*sin(rlat)*cos(rlat)*sin(rlon)
 
      end

**ETIDE0 -- Solid earth tidal displacement caused by one planet.
*+
      FUNCTION ETIDE0 (X, XD, MASS, RE, H2LOVE)
      REAL*8 X(3), XD(3), MASS, RE, H2LOVE, ETIDE0
*
* Compute radial solid earth tidal displacement according to Love model.
*
* Arguments:
*  X       (input): Position of observer in XYZ (earth centered).
*  XD      (input): Position of disturbing body in XYZ (earth centered).
*  MASS    (input): Mass of the disturbing body (Earth = 1).
*  RE      (input): Mean earth radius in meters.
*  H2LOVE  (input): Love number h2.
*  ETIDE0 (output): Radial solid earth tidal displacement in meters.
*-
* 12-Nov-1991: Created (Remko Scharroo).
*-----------------------------------------------------------------------
      real*8 scprod,r,rd,fact

      r=dsqrt(x(1)**2+x(2)**2+x(3)**2)
      rd=dsqrt(xd(1)**2+xd(2)**2+xd(3)**2)
      fact=re*mass*(re/rd)**3*(h2love/2)
      scprod=(x(1)*xd(1)+x(2)*xd(2)+x(3)*xd(3))/r/rd
      etide0=fact*(3*scprod*scprod-1)
      end
