* file: quasob.f
*
* QUASOB is the subroutine that determines the quasi-observation for
* the calculated crossover location by first swichting from longitude
* to time for the x axis
*
      subroutine quasob(jtrack,rdat,rlatx,utcx,seahx,intorb,orbhx)
      integer*4 jtrack,k
      real*8 rlatc,rlatx,rlatxc,ux,utcx,seahx,orbhx,arglat,f1
      logical intorb
      include "maxcom.inc"
      integer*4 mpoint,mplus1,mwork
      parameter (mpoint=6,mplus1=4)
      parameter (mwork=6*mplus1+2*mpoint)
      real*8 x(mpoint),y(mpoint),u(mpoint),rdat(5,mpoint),
     .ss(mplus1),cf(mplus1),work(mwork),
     .aa(mpoint),v1(mpoint),vn(mpoint),cc(mpoint),d(mpoint),
     .incl,orb(mpoint)
*
* Store in temporary arrays: x=time, y=sea surface height, orb=orb height
* (all relative to the third point)
*
      do k=1,mpoint
         x(k)=rdat(1,k)-rdat(1,3)
	 orb(k)=rdat(4,k)-rdat(4,3)
         y(k)=rdat(5,k)-rdat(5,3)
      enddo
*
* Determine epoch UTCX of crossover for current track (JTRACK), by
* assuming a linear relation between time UTC and argument of latitude U.
*
      incl=inc(jtrack)
      do k=2,5
         call geocen(rdat(2,k),rdat(4,k)*1d3,rlatc,f1)
         u(k)=arglat(rlatc,incl)
      enddo
      call e02adf(x(2),u(2),4,2,ss,cf,work)
      call geocen(rlatx,(rdat(4,3)+rdat(4,4))*5d2,rlatxc,f1)
      ux=arglat(rlatxc,incl)
      utcx=(ux-cf(1))/cf(2)
*
* Start calculation of pseudo-observation using a quadratic polynomial
* or cubic spline interpolation for the seaheight for current track
*
      if (specif(2:2).eq.'C') then
         call e02baf(mpoint,-1,y,aa,x,utcx,seahx,f1,v1,vn,cc,d)
         call e02baf(mpoint,-2,y,aa,x,utcx,seahx,f1,v1,vn,cc,d)
	 if (intorb) then
            call e02baf(mpoint,-1,orb,aa,x,utcx,orbhx,f1,v1,vn,cc,d)
            call e02baf(mpoint,-2,orb,aa,x,utcx,orbhx,f1,v1,vn,cc,d)
	 endif
      else if (specif(2:2).eq.'P') then
         call e02adf(x,y,mpoint,4,ss,cf,work)
         seahx=cf(1)+cf(2)*utcx+cf(3)*utcx**2+cf(4)*utcx**3
	 if (intorb) then
            call e02adf(x,orb,mpoint,4,ss,cf,work)
            orbhx=cf(1)+cf(2)*utcx+cf(3)*utcx**2+cf(4)*utcx**3
	 endif
      else
         call e02adf(x,y,mpoint,3,ss,cf,work)
         seahx=cf(1)+cf(2)*utcx+cf(3)*utcx**2
	 if (intorb) then
            call e02adf(x,orb,mpoint,3,ss,cf,work)
            orbhx=cf(1)+cf(2)*utcx+cf(3)*utcx**2
	 endif
      endif
      utcx=utcx+rdat(1,3)
      orbhx=orbhx+rdat(4,3)
      seahx=seahx+rdat(5,3)
      end
