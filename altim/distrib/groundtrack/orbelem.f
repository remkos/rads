      program orbelem

************************************************************************
* Compute orbital elements from ODR files.
*
* Last updates:
* 27-Feb-1997 - Implement xODR
************************************************************************

      implicit none
      integer mrec
      parameter (mrec=500000)
      real*8 year,secday,pi,repeat,rad,ogdot,wmdot,
     |	t(mrec),r(mrec),xyz(3,mrec),lon(mrec),lat(mrec),orb(mrec)
      integer*4 irep,iarc,nrec,remid,ncycle,i
      real*8 og0(mrec),incl(mrec),wm0(mrec),trel(mrec),
     |  incl1,wm1,wmdot1,og1,ogdot1,dum1,dum2
      real*8 v(3),a,vnorm,ort1(3),ort2(3),x,y,
     |  wm,scaprd,w,e,og,anmtot,sfdist
      real*8 mat(6)/6*0d0/,res(3)/3*0d0/,d
      real*8 coswm,sinwm,rlon,rlat,rorb,rr,p,theta
      real*8 s,srms/0d0/,smin/1d40/,smax/-1d40/,rev
      integer period,start,end,fd,tstep,odrinfo,getorb
      integer h(3)/0,1,3/,kp(3),if,n,m,recstep/5/
      character satnm*8,filenm*80,arg*80

      parameter (year=365.24219d0,secday=86400d0)
      parameter (ncycle=501,repeat=35d0)
*
* Set default rates
*
      pi=4*atan(1d0)
      rad=pi/180
      period=nint(5.5d0*86400)
      ogdot=-2*pi/secday
      wmdot=ncycle*2*pi/(repeat*secday)

      call getarg(1,filenm)
      call getarg(2,arg)
      if (arg.eq.' ') then
         i=odrinfo(fd,filenm,satnm,irep,iarc,remid,nrec,
     |		start,end,tstep,x,y,rev)
         call closef(fd)
	 if (i.lt.0) goto 9999
      else
	 read (arg,*) start
	 call getarg(3,arg)
	 read (arg,*) end
	 start=(start-46066)*86400
	 end=(end-46066)*86400
	 nrec=(end-start)/60+1
      endif
      do n=1,nrec
	 t(n)=start+(n-1)*60
	 i=getorb(t(n),lat(n),lon(n),orb(n),'+'//filenm,.true.)
	 if (i.gt.0) goto 9999
      enddo
*
* Convert to Cartesian
*
      do n=1,nrec
	 og=(t(n)-t(1))*ogdot
         call geoxyz(lat(n)*rad,lon(n)*rad-og,orb(n),xyz(1,n),r(n))
      enddo

      m=0
      do n=1+recstep,nrec-recstep,recstep
	 m=m+1
	 trel(m)=t(n)-t(1)
*
* Compute orbit normal
*
         call vecprd(xyz(1,n-recstep),xyz(1,n+recstep),v)
	 a=vnorm(v)
	 do i=1,3
	    v(i)=v(i)/a
	 enddo
	 incl(m)=acos(v(3))
	 og0(m)=atan2(v(1),-v(2))
*
* Compute ortogonal in-orbit vectors
*
	 ort1(1)=cos(og0(m))
	 ort1(2)=sin(og0(m))
	 ort1(3)=0d0
	 ort2(1)=-cos(incl(m))*sin(og0(m))
	 ort2(2)= cos(incl(m))*cos(og0(m))
	 ort2(3)= sin(incl(m))
*
* Make projections along these vectors
*
	 x=scaprd(xyz(1,n),ort1)
	 y=scaprd(xyz(1,n),ort2)
	 wm=atan2(y,x)
	 wm0(m)=wm-wmdot*(t(n)-t(1))
	 i=nint(wm0(m)/2/pi)
	 wm0(m)=wm0(m)-i*2*pi

	 coswm=cos(wm)
	 sinwm=sin(wm)
	 mat(1)=mat(1)+1d0
	 mat(2)=mat(2)+coswm
	 mat(3)=mat(3)+coswm*coswm
	 mat(4)=mat(4)+sinwm
	 mat(5)=mat(5)+sinwm*coswm
	 mat(6)=mat(6)+sinwm*sinwm

	 d=r(n)
         res(1)=res(1)+d
	 res(2)=res(2)+d*coswm
	 res(3)=res(3)+d*sinwm

*	 write (6,*) xyz(1,n),xyz(2,n),xyz(3,n),r(n)
*	 write (6,'(6f12.6)') v,incl(m)/rad,og0(m)/rad,wm0(m)/rad
      enddo
*
* Do linear regression on og0, and wm0
*
      call regres(m,trel,og0,og1,ogdot1,dum1)
*     write (6,*) og1/rad,ogdot1/rad*secday,dum1
      call regres(m,trel,wm0,wm1,wmdot1,dum1)
*     write (6,*) wm1/rad,wmdot1/rad*secday,dum1
      wmdot=wmdot+wmdot1
      ogdot=ogdot+ogdot1
*
* Do statistics on incl
*
      call statis(m,incl,incl1,dum1,dum2)
*     write (6,*) incl1/rad,dum2/rad
*
* Compute a, e, and w
*
      call lincho(mat,h,res,v,kp,3,dum1,if)
      a=res(1)
      e=sqrt(res(2)**2+res(3)**2)/res(1)
      w=atan2(-res(3),-res(2))

      p=a*(1-e*e)
      m=0
      do n=1,nrec,recstep
	 m=m+1
	 wm=wm1+wmdot*(t(n)-t(1))
	 theta=anmtot(wm-w,e)
	 og=og1+ogdot*(t(n)-t(1))
         v(1)=p/(1+e*cos(theta))
	 v(2)=0
	 v(3)=0
	 call rotate(3,-(w+theta),v,v)
	 call rotate(1,-incl1,v,v)
	 call rotate(3,-og,v,v)
	 call xyzgeo(v,rr,rlat,rlon,rorb)
         s=sfdist(rlat,rlon,lat(n)*rad,lon(n)*rad)*6378.137d0
	 smin=min(smin,s)
	 smax=max(smax,s)
	 srms=srms+s*s
	 rlat=rlat/rad
	 rlon=rlon/rad
	 if (rlon.lt.0) rlon=rlon+360
*	 write (6,600) t(n)/86400+46066,rlat,rlon,rorb
*	 write (6,610) rlat,rlon,rorb,lat(n),lon(n),orb(n)
*	 write (6,610) rlat-lat(n),rlon-lon(n),rorb-orb(n),s
      enddo
  600 format (3f12.6,f12.3)
  610 format (2f12.6,f12.3,2f12.6,f12.3)

      write (6,590) iarc,
     |    t(1)/86400+46066,t(nrec)/86400+46066,recstep*60d0,
     |    a,e,incl1/rad,og1/rad,
     |	  ogdot/rad*secday,w/rad,0d0,wm1/rad,
     |	  wmdot/rad*secday,smax,sqrt(srms/m)
  590 format (i3.3,2f9.2,f8.1,f15.6,f15.12,7f15.9,f7.3,f8.3)
9999  end
