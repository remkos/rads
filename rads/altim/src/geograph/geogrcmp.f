*+GEOGRCMP - Compute geographically correlated radial orbit error per component
**
      SUBROUTINE GEOGRCMP (LAT, LON, DRLM_C, DRLM_S)
      REAL*8     LAT, LON, DRLM_C(*), DRLM_S(*)
*
* Compute each component of the geographically correlated mean and
* variable part of radial orbit error, as defined by Rosborough
* (5.42-43). This means, compute the contribution of each C or S
* coefficient of the gravity field.
*
* Arguments:
*  LAT        (input): Geocentric latitude (rad)
*  LON        (input): Longitude (rad)
*  DRLM_C(N) (output): Mean geographically correlated orbit error
*                      for coefficient N
*  DRLM_S(N) (output): Geographically anti-correlated orbit error
*                      for coefficient N
*-
* 22-Feb-1995 - Created from GEOGRAPH
*  9-Aug-2001 - Sped up computation of sin(m*lon) and cos(m*lon)
*----------------------------------------------------------------------
      include "geograph.inc"
      integer l,m,i,ip
      real*8 lat0/-1d40/,lon0/-1d40/
      save lat0,lon0

* For new latitude: initialize Y, Psi, Phi, and Q

      if (lat.ne.lat0) then
         call q_lm(lat)
         lat0=lat
      endif

* For new longitude: initialize cos/sin of multiples of lon

      if (lon.ne.lon0) then
         call mangle(lmax,lon,cosml,sinml)
         lon0=lon
      endif

      do ip=1,ipmax
         l=ideg(ip)
         m=iord(ip)
         i=h(l)+m
	 if (ics(ip).eq.1) then
            drlm_c(ip)=+qc(i)*cs(ip)*cosml(m)
            drlm_s(ip)=+qs(i)*cs(ip)*sinml(m)
         else
            drlm_c(ip)=+qc(i)*cs(ip)*sinml(m)
            drlm_s(ip)=-qs(i)*cs(ip)*cosml(m)
	 endif
      enddo
      end
