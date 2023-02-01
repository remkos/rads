**EPHARG -- Determine ephemeris arguments of sun and moon.
*+
      FUNCTION EPHARG (K, TIME)
      INTEGER*4 K
      REAL*8  TIME, EPHARG
*
* Compute mean longitude of Moon and Sun, their perigee and node,
* or their rates.
* Based on Newcomb' and Brown' Theory.
*
* Angles describing long-term motion of the Sun and the Moon
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
*  K       (input): Argument identifier:
*                   K=1: Greenwich siderial time (theta-g)
*                   K=2: Mean longitude of the Moon (s)
*                   K=3: Mean longitude of the Sun (h)
*                   K=4: Mean longitude of the lunar perigee (p)
*                   K=5: Mean longitude of the lunar node (N)
*                   K=6: Mean longitude of the solar perigee (ps)
*                   K<0: Rate of the argument identified by -K
*  TIME    (input): MJD of the epoch at which argument or argument rate must
*                   be determined.
*  EPHARG (output): Ephemeris argument identified by K (in degrees) or
*                   the rate of the argument identified by -K (in deg/day)
*-
* 12-Nov-1991: Created (Remko Scharroo)
*-----------------------------------------------------------------------
      real*8 ephar0
      if (k.eq.1) then
        epharg=dmod(ephar0(3,time)+dmod(time-.5d0,1d0)*360d0,360d0)
      else if (k.gt.0) then
	 epharg=dmod(ephar0(k,time),360d0)
      else if (k.eq.-1) then
	 epharg=ephar0(-3,time)+360d0
      else
	 epharg=ephar0(k,time)
      endif
      end

      function ephar0(k,time)
      integer*4 k
      real*8    ephar0,time
      real*8 arg(0:3,2:6),tref,t
      parameter (tref=15019.5d0)

      data arg / 270.43659d0,+481267.89057d0,+0.00198d0,+0.000002d0,
     |		 279.69660d0, +36000.76892d0,+0.00030d0, 0.000000d0,
     |		 334.32956d0,  +4069.03403d0,-0.01032d0,-0.000010d0,
     |		 259.18328d0,  -1934.14201d0,+0.00208d0,+0.000002d0,
     |		 281.22083d0,     +1.71902d0,+0.00045d0,+0.000003d0 /
      save arg

      t=(time-tref)/36525
      if (k.gt.0) then
	 ephar0=arg(0,k)+arg(1,k)*t+arg(2,k)*t*t+arg(3,k)*t*t*t
      else
	 ephar0=(arg(1,-k)+arg(2,-k)*2*t+arg(3,-k)*3*t*t)*36525
      endif
      end
