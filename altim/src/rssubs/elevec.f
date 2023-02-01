**ELEVEC -- Convert Keplerian elements to position-velocity vector
*+
      SUBROUTINE ELEVEC (ELEM, VECT, GM)
      REAL*8 ELEM(6), VECT(*), GM
*
* This routine converts the 6 Keplerian elements (a, e, i, O, w, th) to
* a position vector (X, Y, Z) or a position-velocity vector
* (X, Y, Z, Xdot, Ydot, Zdot). When the given gravitational parameter GM is 
* less or equal to zero, then only the position vector is determined.
*
* Arguments:
*  ELEM  (input): Array of 6 Keplerian elements:
*                 ELEM(1): Semi-major axis (meters),
*                 ELEM(2): Eccentricity (unity),
*                 ELEM(3): Inclination (radians),
*                 ELEM(4): Right ascension of the ascending node (radians),
*                 ELEM(5): Argument of the pericenter (radians),
*                 ELEM(6): True anomaly (radians).
*  VECT (output): Array of position (and velocity):
*                 VECT(1)..VECT(3): Position vector (X,Y,Z) (meters),
*                 VECT(4)..VECT(6): Velocity vector (Xdot,Ydot,Zdot)
*                                   (meter/sec).
*  GM    (input): Gravitational parameter of the central body
*                 (meter**3/sec**2). If GM<=0 then no velocity vector is
*                 computed.
*-
* 12-Nov-1991: Created (Remko Scharroo)
*  4-Aug-2004: Allow overriding of input values
*-----------------------------------------------------------------------
      real*8 e,p,cv,ecv,r,u,cu,su,ci,si,cocu,sosu,socu,cosu,
     |	fx,fy,fz,ff,vr,vu,so,co

      e=elem(2)
      p=elem(1)*(1-e**2)
      cv=cos(elem(6))
      ecv=1+e*cv
      r=p/ecv
      u=elem(5)+elem(6)
      cu=cos(u)
      su=sin(u)
      co=cos(elem(4))
      so=sin(elem(4))
      ci=cos(elem(3))
      si=sin(elem(3))
      cocu=co*cu
      sosu=so*su
      socu=so*cu
      cosu=co*su
      fx=cocu-sosu*ci
      fy=socu+cosu*ci
      fz=su*si
      vect(1)=r*fx
      vect(2)=r*fy
      vect(3)=r*fz
      if (gm.le.0) return
      ff=sqrt(gm/p)
      vr=ff*e*sin(elem(6))
      vu=ff*ecv
      vect(4)=vr*fx-vu*(cosu+socu*ci)
      vect(5)=vr*fy-vu*(sosu-cocu*ci)
      vect(6)=vr*fz+vu*cu*si
      end
