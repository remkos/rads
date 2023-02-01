      subroutine cmdprofi
      include "pim.inc"
      include "cmdprofi.inc"
      real trk0,trk1,t0,t1
      real*8 sec85
      integer it0,it1
      logical pimopt,l

* Call PROFILE if requested

      call pgsave
  335 if (pop1cmd('PROFI',argum)) then
	 ci_track=-2
	 ci_profile=1
	 ls=1
	 lw=1
	 factor=1e-2
	 offset=0
         trk0=-1e30
         trk1=1e30
         t0=minint4
         t1=maxint4
	 angle=1e30
	 direction=0
	 name='mean.xgf'
	 l=pimopt('ci=',argum,ci_profile,ci_track,dum,dum)
	 if (ci_track.eq.-2) ci_track=ci_profile
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
	 l=pimopt('fact=',argum,factor,dum,dum,dum)
	 l=pimopt('off=',argum,offset,dum,dum,dum)
	 l=pimopt('angle=',argum,angle,dum,dum,dum)
         if (pimopt('trk=',argum,trk0,trk1,dum,dum)) then
            if (trk1.ge.1e20) trk1=trk0
         endif
         l=pimopt('t=',argum,t0,t1,dum,dum)
	 if (pimopt('-asc',argum,dum,dum,dum,dum)) direction=+1
	 if (pimopt('-des',argum,dum,dum,dum,dum)) direction=-1
         it0=nint(sec85(0,dble(t0)))
         it1=nint(sec85(0,dble(t1)))
         if (it0.ne.minint4 .and. it1.eq.maxint4) it1=it0+86400
	 call strip(argum,name)
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call profile(name,
     .		nint(trk0),nint(trk1),it0,it1,maxgrd/3,work1)
	 goto 335
      endif
      call pgunsa
      end

**PROFILE -- Plot XGF data set.
*+
      SUBROUTINE PROFILE (FILENM,
     .		TRK0, TRK1, T0, T1, NWORK, WORK)
      CHARACTER*(*) FILENM
      INTEGER NWORK, TRK0, TRK1, T0, T1
      REAL WORK(NWORK,3)
*
* Routine plots an XGF data set FILENM as a profile in a plot. It buffers
* a number of NWORK points in the buffer WORK first.
*
* Arguments:
*  FILENM  (input): XGF filenm
*  TRK0, TRK1 (in): Limits on tracknumber.
*  T0, T1  (input): Limits on time (Sec85).
*  NWORK   (input): Size of the workspace WORK.
*  WORK           : Working space.
*-
      integer i,n,lat,lat0,lon,lon0,nrec,itime,itime0,unit,hgt,maxint4
      integer*2 sig,trk
      integer dum,freeunit
      character spec*4,spec2*2
      parameter (maxint4=2000000000)
*
* Open XGF file on new unit
*
      write (6,550) 'Reading and plotting XGF profile ...'
550   format (a)
      unit=freeunit()
      open (unit,file=filenm,status='old',form='unformatted',
     .   access='direct',recl=18)
      read (unit,rec=1) spec,nrec,dum,dum,spec2
      if (spec.ne.'@XGF') return

      itime0=0
      n=0
      do i=1,nrec
	 read (unit,rec=i+1,err=50) itime,lat,lon,hgt,sig
         if (spec2.ne.'Tm') then
            if (itime.lt.t0 .or. itime.gt.t1) goto 50
         else if (sig.lt.0) then
            if (itime.lt.t0 .or. itime.gt.t1) goto 50
         else
            if (trk.lt.trk0 .or. trk.gt.trk1) goto 50
         endif
	 if (n.eq.nwork .or. itime-itime0.gt.10
     .     .or. abs(lon-lon0).gt.500000 .or. abs(lat-lat0).gt.500000)
     .		call profile1(n,work(1,1),work(1,2),work(1,3))
	 n=n+1
	 work(n,1)=lon/1e6
	 work(n,2)=lat/1e6
	 work(n,3)=hgt/1e6
	 itime0=itime
	 lon0=lon
	 lat0=lat
50       continue
      enddo
      call profile1(n,work(1,1),work(1,2),work(1,3))

      end

      subroutine profile1(n,x,y,z)
      integer n,i
      real x(*),y(*),z(*)
      real dx,dy,a,pi
      real xw0,xw1,yw0,yw1,xp0,xp1,yp0,yp1,fx,fy

      include "cmdprofi.inc"

      if (n.le.1) goto 100
      if (direction.eq.-1 .and. y(1).lt.y(n)) goto 100
      if (direction.eq.+1 .and. y(1).gt.y(n)) goto 100

      pi=4*atan(1d0)
      call pmconv(n,x,y)
      if (ci_track.ge.0) then
	 call pgsci(nint(ci_track))
         call pgline(n,x,y)
      endif

      call pgqwin(xw0,xw1,yw0,yw1)
      call pgqvp(2,xp0,xp1,yp0,yp1)
      fx=(xw1-xw0)/(xp1-xp0)*1d3*factor
      fy=(yw1-yw0)/(yp1-yp0)*1d3*factor

      dx=x(2)-x(1)
      dy=y(2)-y(1)

      do i=1,n
	 if (angle.gt.1e10) then
	    a=atan2(dy,dx)
	    if (a.lt.0) a=a+pi
	 else
	    a=-angle*pi/180	! angle is measured from due north
	 endif
	 dx=x(i+1)-x(i)
	 dy=y(i+1)-y(i)
	 x(i)=x(i)-(z(i)-offset)*sin(a)*fx
	 y(i)=y(i)+(z(i)-offset)*cos(a)*fy
      enddo
      if (ci_profile.ge.0) then
	 call pgsci(nint(ci_profile))
         call pgline(n,x,y)
      endif

100   n=0

      end
