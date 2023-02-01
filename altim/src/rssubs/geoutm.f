**GEOUTM -- Convert geodetic latitude and longitude to UTM coordinates
*+
      SUBROUTINE GEOUTM (CENMED, LAT, LON, X, Y)
      REAL*8 CENMED, LAT, LON, X, Y
*
* This routine converts the geodetic coordinates of a point on the
* Earth to the Universal Transvers Mercator coordinates (UTM) on the
* GRS80 spheroid.
*
* Arguments:
*  CENMED (input) : Central meridian of the projection (radians)
*  LAT    (input) : Geodetic latitude of the point (radians)
*  LON    (input) : East longitude (radians).
*  X     (output) : UTM X-coordinate (meters).
*  Y     (output) : UTM Y-coordiante (meters).
*-
* 17-Jan-1992 - Created [Remko Scharroo]
*-----------------------------------------------------------------------
      INCLUDE 'GRS80.inc'
      real*8 ep2,e2,an,an2,an3,an4,an5,a,b,c,d,e,zkphi,x0,y0,t2,t4,
     |       s,c4,c5,t1,ve,sv,sa,sb,sr,st,su,df,delta,c1,c2,c3,alat,
     |       s1,s2

      parameter (
     |	ep2=(ae**2-ap**2)/ap**2,
     |	e2=(ae**2-ap**2)/ae**2,
     |	an=(ae-ap)/(ae+ap),
     |	an2=an*an,an3=an2*an,an4=an2*an2,an5=an2*an3,
     |	a=ae*(1-an+5/4d0*(an2-an3)+81/64d0*(an4-an5)),
     |	b=3/2d0*ae*(an-an2+7/8d0*(an3-an4)+55/64d0*an5),
     |	c=15/16d0*ae*(an2-an3+3/4d0*(an4-an5)),
     |	d=35/48d0*ae*(an3-an4+11/16d0*an5),
     |	e=315/512d0*ae*(an4-an5),
     |	zkphi=.9996d0,x0=500d3,y0=0d3)

      df=cenmed-lon
      delta=abs(df)
      alat=abs(lat)
      s1=sin(alat)
      s2=s1*s1
      c1=cos(alat)
      c2=c1*c1
      c3=c2*c1
      c4=c2*c2
      c5=c2*c3
      t1=tan(alat)
      t2=t1*t1
      t4=t2*t2

      s=zkphi*(a*alat-b*sin(2*alat)+c*sin(4*alat)-d*sin(6*alat)+
     |  e*sin(8*alat))
      ve=ae/sqrt(1-e2*s2)
      sr=ve/2*s1*c1*zkphi
      st=ve*s1*c3*(5-t2+9*ep2*c2+4*ep2*ep2*c4)*zkphi/24
      su=ve*c1*zkphi
      sv=ve/6*c3*(1-t2+ep2*c2)*zkphi
      sa=delta**6/720*ve*s1*c5*(61-58*t2+
     |	t4+270*ep2*c2-330*ep2*s2)*zkphi
      sb=delta**5/120*ve*c5*(5-18*t2*t4+14*ep2*c2-58*ep2*s2)*zkphi

      x=su*delta+sv*delta**3+sb
      if (df.gt.0) x=-x
      x=x+x0

      y=s+sr*delta**2+st*delta**4+sa
      if (lat.lt.0) y=-y
      y=y+y0

      end
