**VECELE -- Convert state vector to Keplerian elements
*+
      SUBROUTINE VECELE (VECT, ELEM, GM)
      REAL*8 VECT(6), ELEM(6), GM
*
* This routine converts a state-vector (X, Y, Z, Xdot, Ydot, Zdot) to
* the 6 Keplerian elements (a, e, i, O, w, th).
* The gravitational parameter GM has to be given as well.
*
* Arguments:
*  VECT  (input): Array of position (and velocity):
*                 VECT(1)..VECT(3): Position vector (X,Y,Z) (meters),
*                 VECT(4)..VECT(6): Velocity vector (Xdot,Ydot,Zdot)
*                                   (meter/sec).
*  ELEM (output): Array of 6 Keplerian elements:
*                 ELEM(1): Semi-major axis (meters),
*                 ELEM(2): Eccentricity (unity),
*                 ELEM(3): Inclination (radians),
*                 ELEM(4): Right ascension of the ascending node (radians),
*                 ELEM(5): Argument of the pericenter (radians),
*                 ELEM(6): True anomaly (radians).
*  GM    (input): Gravitational parameter of the central body
*                 (meter**3/sec**2). If GM<=0 then no velocity vector is
*                 computed.
*-
*  2-May-1995: Created (Remko Scharroo)
*-----------------------------------------------------------------------
      include "math.inc"

      real*8 a,e,rincl,rnode,omega,true,aeiowf(6)
      equivalence (a,aeiowf(1)),(e,aeiowf(2)),(rincl,aeiowf(3))
      equivalence (rnode,aeiowf(4)),(omega,aeiowf(5)),(true,aeiowf(6))
      real*8 x,y,z,xdot,ydot,zdot,xyzxyz(6)
      equivalence (x,xyzxyz(1)),(y,xyzxyz(2)),(z,xyzxyz(3))
      equivalence (xdot,xyzxyz(4)),(ydot,xyzxyz(5)),(zdot,xyzxyz(6))
      integer*4 i
      real*8 rrdot,r,vsq,ainv
      real*8 hx,hy,hz,hsini2,hsq,h,ome2
      real*8 cosi,sini,sini2,hsini
      real*8 resinf,recosf
      real*8 suprod,cuprod

* Store input vector in internal array xyzxyz (equivalenced to x,y,z,
* xdot,ydot,zdot)

      do i=1,6
	 xyzxyz(i)=vect(i)
      enddo

* Compute semi-major axis

      rrdot=x*xdot+y*ydot+z*zdot
      r=sqrt(x**2+y**2+z**2)
      vsq=xdot**2+ydot**2+zdot**2
      ainv=two/r-vsq/gm
      a=one/ainv

* Compute eccentricity

      hx=y*zdot-z*ydot
      hy=z*xdot-x*zdot
      hz=x*ydot-y*xdot
      hsini2=hx**2+hy**2
      hsq=hsini2+hz**2
      h=sqrt(hsq)
      ome2=hsq*ainv/gm
      e=sqrt(one-ome2)

* Compute inclination

      cosi=hz/h
      sini2=one-cosi**2
      sini=sqrt(sini2)
      hsini=h*sini
      rincl=atan2(sini,cosi)
      if (rincl.lt.zero) rincl=rincl+twopi

* Compute longitude of ascending node

      rnode=atan2(hx,-hy)
      if (rnode.lt.zero) rnode=rnode+twopi

* Compute true anomaly

      resinf=a*ome2*rrdot/h
      recosf=a*ome2-r
      true=atan2(resinf,recosf)
      if (true.lt.zero) true=true+twopi

* Compute argument of perigee

      suprod=-hz*(x*hx+y*hy)+z*hsini2
      cuprod=h*(-x*hy+y*hx)
      omega=atan2(suprod*recosf-cuprod*resinf,
     |		  cuprod*recosf+suprod*resinf)
      if (omega.lt.zero) omega=omega+twopi

* Export aeiowf (equivalenced to a,e,rincl,rnode,omega,true) to ele

      do i=1,6
	 elem(i)=aeiowf(i)
      enddo
      end
