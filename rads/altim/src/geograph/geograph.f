*+GEOGRAPH - Compute geographically correlated radial orbit error
**
      SUBROUTINE GEOGRAPH (LAT, LON, DR_C, DR_S)
      REAL*8     LAT, LON, DR_C, DR_S
*
* Compute geographically correlated mean and variable part of radial orbit
* error, as defined by Rosborough (5.42-43) or the the associated
* variances, as defined by Rosborough (7.33-34).
*
* Arguments:
*  LAT     (input): Geocentric latitude (rad)
*  LON  (input): Longitude (rad)
*  DR_C   (output): Mean geographically correlated orbit error
*  DR_S   (output): Geographically anti-correlated orbit error
*-
* 18-Mar-1993 - Created (Remko Scharroo)
*  2-Feb-1994 - Computation of variances added.
*  9-Aug-2001 - Sped up computation of sin(m*lon) and cos(m*lon)
*----------------------------------------------------------------------
      include "geograph.inc"
      integer i,j,k,l,m,ip,jp
      real*8 lat0/-1d40/,lon0/-1d40/
      real*8 pcc(2,2,0:ndeg,0:ndeg),pss(2,2,0:ndeg,0:ndeg)
      save lat0,lon0,pcc,pss

* For new latitude: initialize Y, Psi, Phi, and Q

      if (lat.ne.lat0) then
         if (test) write (*,*) 'Initialize Y, Psi, Phi and Q'
         call q_lm(lat)
         if (docovar) then
            do k=0,lmax
               do m=0,lmax
                  pcc(1,1,m,k)=0
                  pcc(2,1,m,k)=0
                  pcc(1,2,m,k)=0
                  pcc(2,2,m,k)=0
                  pss(1,1,m,k)=0
                  pss(2,1,m,k)=0
                  pss(1,2,m,k)=0
                  pss(2,2,m,k)=0
               enddo
            enddo
            do ip=1,ipmax
               do jp=1,ip
                  l=ideg(ip)
                  m=iord(ip)
                  j=ideg(jp)
                  k=iord(jp)
                  pcc(ics(ip),ics(jp),m,k)=pcc(ics(ip),ics(jp),m,k)+
     |               qc(h(l)+m)*qc(h(j)+k)*covar(g(ip)+jp)
                  pss(ics(ip),ics(jp),m,k)=pss(ics(ip),ics(jp),m,k)+
     |               qs(h(l)+m)*qs(h(j)+k)*covar(g(ip)+jp)
               enddo
            enddo
         endif
         lat0=lat
      endif

* For new longitude: initialize cos/sin of multiples of lon

      if (lon.ne.lon0) then
         call mangle(lmax,lon,cosml,sinml)
         lon0=lon
      endif

* Initialize mean and variable part of gravity induced radial orbit
* error

      dr_c=0
      dr_s=0

* If covariances are given:
* Note: Because pcs(n) contains only half the complete variance-covariance
* matrix, pcs(n) are pre-multiplied by the appropriate factors, i.e.
* * 2, if CS covariance (to compensate for missing SC)
* * 2, if CC or SS covariance and lm <> jk (to compensate for missing
*                  covariance CC and SS with lm < jk)
* * 1, if CC or SS variance with lm = jk (diagonal is complete) 

      if (docovar) then
         do k=0,lmax
            do m=0,lmax
               dr_c=dr_c+
     |            pcc(1,1,m,k)*cosml(m)*cosml(k)+
     |            pcc(1,2,m,k)*cosml(m)*sinml(k)+
     |            pcc(2,1,m,k)*sinml(m)*cosml(k)+
     |            pcc(2,2,m,k)*sinml(m)*sinml(k)
               dr_s=dr_s+
     |            pss(1,1,m,k)*sinml(m)*sinml(k)-
     |            pss(1,2,m,k)*sinml(m)*cosml(k)-
     |            pss(2,1,m,k)*cosml(m)*sinml(k)+
     |            pss(2,2,m,k)*cosml(m)*cosml(k)
            enddo
         enddo
         dr_c=sqrt(dr_c)
         dr_s=sqrt(dr_s)

* If coefficients are given

      else
         do ip=1,ipmax
            l=ideg(ip)
            m=iord(ip)
            i=h(l)+m
            if (ics(ip).eq.1) then
               dr_c=dr_c+qc(i)*cs(ip)*cosml(m)
               dr_s=dr_s+qs(i)*cs(ip)*sinml(m)
            else
               dr_c=dr_c+qc(i)*cs(ip)*sinml(m)
               dr_s=dr_s-qs(i)*cs(ip)*cosml(m)
            endif
         enddo
      endif
      end
