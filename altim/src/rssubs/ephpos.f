**EPHPOS -- Determine ephemeris of Sun and Moon.
*+
      SUBROUTINE EPHPOS (PLANET, TIME, XYZ)
      INTEGER*4 PLANET
      REAL*8  TIME, XYZ(3)
*
* Compute the inertial position of the Moon and Sun.
* Based on Newcomb' and Brown' Theory.
*
* The XYZ coordinates are inferred from the angles describing the long-term
* motion of the Sun and the Moon
*  a = a0 + a1*T + a2*T**2 + a3*T**3  with
*  T = centuries past JD 2415020 (31.5 Dec 1899) = 15019.5 MJD
*
*    a0            a1          a2         a3
* 270.43659   +481267.89057  +0.00198  +0.000002  s = mean long. of the Moon
* 279.69660    +36000.76892  +0.00030   0.000000  h = mean long. of the Sun
* 334.32956     +4069.03403  -0.01032  -0.000010  p = mean long. of lunar perig
* 259.18328     -1934.14201  +0.00208  +0.000002  N = mean long. of lunar node
* 281.22083        +1.71902  +0.00045  +0.000003  ps= mean long. of solar perig
*
* Ref: Stefano Casotto - Nominal ocean tide models for TOPEX precise orbit
*      determination, University of Texas at Austin, Dec 1989.
*
* Arguments:
*  PLANET (input): Planet identifier: 1=Moon, 2=Sun.
*  TIME   (input): MJD of the epoch at which argument or argument rate must
*                  be determined.
*  XYZ   (output): X, Y, and Z component of the ephemeris of the Moon or Sun
*                  in meters.
*-
*  9-Feb-1994: Created (Remko Scharroo)
*-----------------------------------------------------------------------
      include 'etide.inc'
      include 'math.inc'
      
      real*8 elem(6),a,anmtot,epharg

      if (planet.eq.1) then
         elem(1)=amoon
         elem(2)=emoon
         elem(3)=imoon*rad
         elem(4)=epharg(5,time)*rad
         a=epharg(4,time)*rad
         elem(5)=a-elem(4)
         elem(6)=anmtot(epharg(2,time)*rad-a,elem(2))
         call elevec(elem(1),xyz,-1d0)
         call rotate(1,-elem(3),xyz,xyz)
      else if (planet.eq.2) then
         elem(1)=asun
         elem(2)=esun
         elem(3)=isun*rad
         elem(4)=0d0
         elem(5)=epharg(6,time)*rad
         elem(6)=anmtot(epharg(3,time)*rad-elem(5),elem(2))
         call elevec(elem(1),xyz,-1d0)
      endif
 
      end
