#ifdef DIFF
*+F_LMPDD -- Initialize normalized inclination functions and derivatives
*
      SUBROUTINE F_LMPDD (LMAX, INCL, PNT, FLMP, DFLMP, DDFLMP)
      INTEGER LMAX, PNT(0:LMAX)
      REAL*8  INCL, FLMP(0:*), DFLMP(0:*), DDFLMP(0:*)
#else
*+F_LMP -- Initialize normalized inclination functions
*
      SUBROUTINE F_LMP (LMAX, INCL, PNT, FLMP)
      INTEGER LMAX, PNT(0:LMAX)
      REAL*8  INCL, FLMP(0:*)
#endif
*
* This subroutine computes normalized inclination functions
* by means of a Fourier analysis on normalized Legendre polynomials.
#ifdef DIFF
* Its first and second derivatives with repect to the inclinations
* are computed as well.
#endif
*
* Input arguments:
*   LMAX      : maximum degree. LMAX may not be greater than 70.
*   INCL      : orbit inclination (radians)
*
* Output arguments:
*   PNT(j)    : Pointer used for FLMP for j=1..LMAX
*   FLMP(j)   : Array of F_lmp         with j=pnt(l)+m*(l+1)+p
#ifdef DIFF
*   DFLMP(j)  : Array of dF_lmp/di     with j=pnt(l)+m*(l+1)+p
*   DDFLMP(j) : Array of d^2F_lmp/di^2 with j=pnt(l)+m*(l+1)+p
#endif
* The matrix FLMP should be dimensioned at least (0:jmax)
* jmax=(2*LMAX+3)*(LMAX+2)*(LMAX+1))/6-1
*-
*  8-Aug-2001 - Created from FINC by Remko Scharroo
*-----------------------------------------------------------------------
      integer ndeg
      parameter (ndeg=70)
      real*8 plm(0:ndeg)
      real*8 cosml(0:ndeg),sinml(0:ndeg),cospu(0:ndeg),sinpu(0:ndeg)
      real*8 cosi,sini,lon,fac,pi,u,du,phi,c,s,cms,cps
      integer l,m,p,i,i1,i2,i2p,par(0:ndeg)
      parameter (pi=3.14159265358979d0)
#ifdef DIFF
      real*8 dc,ds,ddc,dds,cosl,sinl,cosp,sinp,f1,f2,sinl2,sinp2
      real*8 dplm(0:ndeg),ddplm(0:ndeg)
#endif

* Test LMAX against the precompiled maximum NDEG.

      if (lmax.gt.ndeg) stop "f_lmp: lmax > ndeg"

* Make some initializations.

      sini=dsin(incl)
      cosi=dcos(incl)

      do i=0,(2*lmax+3)*(lmax+2)*(lmax+1)/6-1
         flmp(i)=0d0
#ifdef DIFF
         dflmp(i)=0d0
         ddflmp(i)=0d0
#endif
      enddo
      flmp(0)=1d0
      pnt(0)=0
      par(0)=1

* Repeat the following process for each degree L until LMAX.

      do l=1,lmax

* Store the next element of pointer PNT.
* Store the parity of L: PAR(L) is 1 for L is even, -1 for L is odd.
* Determine the stepsize of U.

        pnt(l)=pnt(l-1)+l**2
	par(l)=-par(l-1)
        fac=0.5d0/(l/2+1)
        du=pi*fac

* Cycle through several values of U, starting at dU/2 and increasing
* with dU.

        u=du/2
        do i=0,l/2
          call mangle(l,u,cospu,sinpu)
          phi=dasin(sini*sinpu(1))
          lon=dasin(cosi*sinpu(1)/dcos(phi))
          call mangle(l,lon,cosml,sinml)

* Request the normalized Langrange polynomials for degree L.

#ifdef DIFF
          call p_lmdd(l,phi,plm,dplm,ddplm)

          cosl=cosml(1)
          sinl=sinml(1)
          cosp=dcos(phi)
          sinp=dsin(phi)
	  f1=cosl*sinp/cosp
	  f2=cosl/cosp/cosp
	  sinl2=sinl**2
	  sinp2=sinp**2
#else
          call p_lm(l,phi,plm)
#endif

* Add contribution to F_lmp for each M between 0 and L and for all L
* with (L-P) is even.

          do m=0,l
            i1=par(l-m)
            i2=pnt(l)+m*(l+1)
            c=   plm(m)*fac*(cosml(m)+sinml(m))
            s=i1*plm(m)*fac*(cosml(m)-sinml(m))
#ifdef DIFF
            cps=cosml(m)+sinml(m)
            cms=cosml(m)-sinml(m)
            dc=   fac*(dplm(m)*cps*sinl-plm(m)*m*f1*cms)
            ds=i1*fac*(dplm(m)*cms*sinl+plm(m)*m*f1*cps)
            ddc=   fac*(ddplm(m)*cps*sinl2
     |       +dplm(m)*(-cps*cosl-2*m*cms*sinl)*f1
     |       +plm(m)*(-cms*sinl*(sinp2+1)-m*cps*cosl*sinp2)*m*f2)
            dds=i1*fac*(ddplm(m)*cms*sinl2
     |       +dplm(m)*(-cms*cosl+2*m*cps*sinl)*f1
     |       +plm(m)*(cps*sinl*(sinp2+1)-m*cms*cosl*sinp2)*m*f2)
#endif
            do p=l,0,-2
              i2p=i2+(l+p)/2
              cms=cospu(p)*i1-sinpu(p)
              cps=cospu(p)*i1+sinpu(p)
              flmp(i2p)=flmp(i2p)+c*cms+s*cps
#ifdef DIFF
              dflmp(i2p)=dflmp(i2p)+dc*cms+ds*cps
              ddflmp(i2p)=ddflmp(i2p)+ddc*cms+dds*cps
#endif
* Do not repeat the summation when P is 0.

              if (p.ne.0) then
                i2p=i2p-p
                flmp(i2p)=flmp(i2p)+c*cps+s*cms
#ifdef DIFF
                dflmp(i2p)=dflmp(i2p)+dc*cps+ds*cms
                ddflmp(i2p)=ddflmp(i2p)+ddc*cps+dds*cms
#endif
              endif
            enddo
          enddo
          u=u+du
        enddo
      enddo

      end
