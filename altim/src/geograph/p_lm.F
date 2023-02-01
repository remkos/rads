#ifdef DIFF
**P_LMDD -- Compute Legendre polynomials of order L
*+
      SUBROUTINE P_LMDD (L, PHI, PLM, DPLM, DDPLM)
      INTEGER L
      REAL*8  PHI, PLM(0:*), DPLM(0:*), DDPLM(0:*)
#else
**P_LM -- Compute Legendre polynomials of order L
*+
      SUBROUTINE P_LM (L, PHI, PLM)
      INTEGER L
      REAL*8  PHI, PLM(0:*)
#endif
*
* Computation of normalized Legendre polynomials P_lm(sin PHI)
* for latitude PHI for a single degree L. L may not be greater than 70.
#ifdef DIFF
* The first and second derivative of P_lm with respect to latitude is
* computed as well.
#endif
*
* Arguments:
*  L    (input): Degree for with Legendre polynomials are to be evaluated
*  PHI  (input): Latitude in radians
*  PLM (output): Array of P_lm for each degree M, M=0...L
#ifdef DIFF
*  PLMD  (output): Array of dP_lm/dphi for each degree M, M=0...L
*  PLMDD (output): Array of d^2P_lm/dphi^2 for each degree M, M=0...L
#endif
*-
*  8-Aug-2001 - Created out of LEGDEG by Remko Scharroo
*-----------------------------------------------------------------------
      integer m,ndeg
      parameter (ndeg=360)
      real*8  sinp,cosp,tanp,root(0:2*ndeg+1)
#ifdef DIFF
      real*8  f1,f2
#endif
      logical first/.true./
      save    root,first

* Test degree L against the precompiled maximum NDEG

      if (l.gt.ndeg) stop "p_lm: l > ndeg"

* Initialize lookup table for sqrt(m) on first call only

      if (first) then
         do m=0,2*ndeg+1
            root(m)=dsqrt(dfloat(m))
         enddo
         first=.false.
      endif

* Store sin(PHI), cos(PHI) and tan(PHI)

      sinp=dsin(phi)
      cosp=dcos(phi)
      tanp=sinp/cosp

* Initialize P_00 and then compute successive P_mm using a recursive
* relation. At last P_ll is found.

      plm(0)=1d0
      plm(1)=root(3)*cosp
#ifdef DIFF
      dplm(0)=0d0
      dplm(1)=-root(3)*sinp
      ddplm(0)=0d0
      ddplm(1)=-root(3)*cosp
#endif
      do m=2,l
#ifdef DIFF
         f1=root(2*m+1)/root(2*m)
         plm(m)=f1*cosp*plm(m-1)
         dplm(m)=f1*(-sinp*plm(m-1)+cosp*dplm(m-1))
         ddplm(m)=f1*(-cosp*plm(m-1)-2*sinp*dplm(m-1)
     |      +cosp*ddplm(m-1))
#else
         plm(m)=root(2*m+1)/root(2*m)*cosp*plm(m-1)
#endif
      enddo

* Compute P_l,l-1 from P_l-1,l-1

#ifdef DIFF
      f1=root(2*l+1)
      ddplm(l-1)=f1*(-sinp*plm(l-1)+2*cosp*dplm(l-1)
     |    +sinp*ddplm(l-1))
      dplm(l-1)=f1*(cosp*plm(l-1)+sinp*dplm(l-1))
      plm(l-1)=f1*sinp*plm(l-1)
#else
      plm(l-1)=root(2*l+1)*sinp*plm(l-1)
#endif

* If L is less than 2 we end here

      if (l.lt.2) return

* Use a recursive relation to determine P_l,m from P_l,m+1 and P_l,m+2.
* Start with the computation of P_l,l-2 and go down to P_l,0.

      do m=l-2,0,-1
#ifdef DIFF
         f1=-m*tanp
         f2=root(l-m)*root(l+m+1)
         plm(m)=(2*(m+1)*tanp*plm(m+1)
     |     -root(l-m-1)*root(l+m+2)*plm(m+2))/f2
         dplm(m)=f1*plm(m)+f2*plm(m+1)
         ddplm(m)=f1*dplm(m)-m/cosp**2*plm(m)+f2*dplm(m+1)
#else
         plm(m)=(2*(m+1)*tanp*plm(m+1)
     |     -root(l-m-1)*root(l+m+2)*plm(m+2))/root(l-m)/root(l+m+1)
#endif
      enddo
      plm(0)=plm(0)/root(2)
#ifdef DIFF
      dplm(0)=dplm(0)/root(2)
      ddplm(0)=ddplm(0)/root(2)
#endif
      end
