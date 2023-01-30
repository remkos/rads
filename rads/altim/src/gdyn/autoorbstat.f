**AUTOORBSTAT -- A program to validate the automatic generated orbits
*
*-
*  8-Jan-1999 - Created by Remko Scharroo
*  8-Feb-1999 - Changed grey bars for manoeuvres.
*  3-Nov-1999 - Changed sec85(2,) to sec85(0,)
*-----------------------------------------------------------------------
      program	lodr
      integer	mstep,mman,nday,nstep,nman
      parameter (mstep=1000*1440,mman=200)
      integer	time0,time1,tstep
      integer	i,j,ios,getorb,t
      integer	yymmdd,mdate,man_date(mman),man_length(mman)
      character arg*80,date0*30,date1*30,ascout*80,string*6
      real*8	lat,lon,orb(mstep,5),sec85,stat(mstep,4)
      real	x0,y0,y1,z0/0.0/,z1/60.0/
      logical	verbose/.false./

* Read arguments

      call getarg(1,date0)
      read (date0,*) lat
      time0=nint(sec85(0,lat))
      call getarg(2,date1)
      read (date1,*) lat
      time1=nint(sec85(0,lat))
      call getarg(3,arg)
      read (arg,*) tstep
      nstep=(time1-time0)/tstep
      if (nstep.gt.mstep) call fin('too many timesteps')
      nday=(time1-time0)/86400
      call getarg(5,ascout)

* Read manoeuvre file

      open (10,file='/user/remko/ers/ODR.ERS-2/maneuver.txt')
      do i=1,14
         read (10,*)
      enddo
      nman=0
10    read (10,'(5x,i6,9x,i5)',end=20) yymmdd,i
      nman=nman+1
      man_date(nman)=mdate(2,yymmdd)-46066
      man_length(nman)=i
      goto 10
20    continue

* Interpolate the orbits, start to finish
* 1) Automated orbit, prediction

      do t=time0,time1,tstep
	 call strf1985(arg,
     |   '+/user/remko/ers/ODR.ERS-2/dgm-e04.auto/ODR.%y%m%d',t-86400)
	 i=(t-time0)/tstep+1
	 ios=getorb(dble(t),lat,lon,orb(i,1),arg,verbose)
	 if (ios.gt.0) orb(i,1)=1d30
      enddo

* 2) Automated orbit, restituted

      do t=time0,time1,tstep
	 call strf1985(arg,
     |   '+/user/remko/ers/ODR.ERS-2/dgm-e04.auto/ODR.%y%m%d',t)
	 i=(t-time0)/tstep+1
	 ios=getorb(dble(t),lat,lon,orb(i,2),arg,verbose)
	 if (ios.gt.0) orb(i,2)=1d30
      enddo

* 3) Fast-delivery orbit

      do t=time0,time1,tstep
	 i=(t-time0)/tstep+1
	 arg='/user/remko/ers/ODR.ERS-2/dgm-e04.fd'
	 ios=getorb(dble(t),lat,lon,orb(i,3),arg,verbose)
	 if (ios.gt.0) orb(i,3)=1d30
      enddo

* 4) Preliminary precise orbit

      do t=time0,time1,tstep
	 i=(t-time0)/tstep+1
	 arg='/user/remko/ers/ODR.ERS-2/dgm-e04.prelim'
	 ios=getorb(dble(t),lat,lon,orb(i,4),arg,verbose)
	 if (ios.gt.0) orb(i,4)=1d30
      enddo

* 5) Final precise orbit

      do t=time0,time1,tstep
	 i=(t-time0)/tstep+1
	 arg='/user/remko/ers/ODR.ERS-2/dgm-e04'
         ios=getorb(dble(t),lat,lon,orb(i,5),arg,verbose)
	 if (ios.gt.0) orb(i,5)=1d30
      enddo

* Open plot device

      call getarg(4,arg)
      call pgbeg(0,arg,1,2)
      call pgscr(3,0.,1.,1.)

*** Start the TOP PLOT

      call pgpage
      call pgsvp(0.03,0.97,0.10,0.97)
      call pgswin(0.0,real(nday),z0,z1/2)

* Plot manoeuvres as grey bars of the right length

      call pgsci(15)
      do i=1,nman
	 x0=(man_date(i)-time0/86400)
	 if (man_date(i).eq.man_date(i-1)) then
	    y0=y1+1
	 else
	    y0=z0
	 endif
	 y1=y0+man_length(i)
         call pgrect(x0,x0+1,y0,y1)
      enddo

* Plot the box

      call pgsci(1)
      call pgbox(' ',0.0,0,'m',0.0,0)
      call pgswin(0.0,real(nday),z0,z1)
      call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
	 call pglab('time from start (days)',
     |		'RMS difference between orbit solutions (cm)',
     |		' ')
      call pgmtxt('R',2.5,0.5,0.5,'length of manoeuvre (s)')
      call pgmtxt('B',3.0,0.0,0.5,date0)
      call pgmtxt('B',3.0,1.0,0.5,date1)
      call pgsch(1.5)
      call pgsci(1)
      call pgmtxt('T',-1.5,.02,0.,'\\(0791) predicted vs restituted')
      call pgsci(2)
      call pgmtxt('T',-3.0,.02,0.,'\\(0791) predicted vs fast-delivery')
      call pgsci(3)
      call pgmtxt('T',-4.5,.02,0.,'\\(0791) predicted vs preliminary')
      call pgsci(4)
      call pgmtxt('T',-6.0,.02,0.,'\\(0791) predicted vs precise')
      call pgsch(1.0)

* Now plot lines for each solution difference
* black = predicted - restituted
* red	= predicted - fast-delivery
* blue  = predicted - preliminary
* green = predicted - precise

      call pgslw(1)
      call pgsci(1)
      call difstat(nstep,tstep,orb(1,1),orb(1,2),stat(1,1))
      call pgsci(2)
      call difstat(nstep,tstep,orb(1,1),orb(1,3),stat(1,2))
      call pgsci(3)
      call difstat(nstep,tstep,orb(1,1),orb(1,4),stat(1,3))
      call pgsci(4)
      call difstat(nstep,tstep,orb(1,1),orb(1,5),stat(1,4))
      call pgslw(1)

*** Start BOTTOM PLOT

      call pgpage
      call pgswin(0.0,real(nday),z0,z1/2)

* Plot manoeuvres as grey bars of the right length

      call pgsci(15)
      do i=1,nman
	 x0=(man_date(i)-time0/86400)
	 if (man_date(i).eq.man_date(i-1)) then
	    y0=y1+1
	 else
	    y0=z0
	 endif
	 y1=y0+man_length(i)
         call pgrect(x0,x0+1,y0,y1)
      enddo

* Plot the box

      call pgsci(1)
      call pgbox(' ',0.0,0,'m',0.0,0)
      call pgswin(0.0,real(nday),z0,z1)
      call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
	 call pglab('time from start (days)',
     |		'RMS difference between orbit solutions (cm)',
     |		' ')
      call pgmtxt('R',2.5,0.5,0.5,'length of manoeuvre (s)')
      call pgmtxt('B',3.0,0.0,0.5,date0)
      call pgmtxt('B',3.0,1.0,0.5,date1)
      call pgsch(1.5)
*     call pgsci(1)
*     call pgmtxt('T',-1.5,.02,0.,'\\(0791) restituted vs predicted')
      call pgsci(2)
      call pgmtxt('T',-3.,.02,0.,'\\(0791) restituted vs fast-delivery')
      call pgsci(3)
      call pgmtxt('T',-4.5,.02,0.,'\\(0791) restituted vs preliminary')
      call pgsci(4)
      call pgmtxt('T',-6.0,.02,0.,'\\(0791) restituted vs precise')
      call pgsch(1.0)

* Now plot lines for each solution difference
* (black = restituted - predicted)
* red	= restituted - fast-delivery
* blue  = restituted - preliminary
* green = restituted - precise

      call pgslw(1)
*     call pgsci(1)
*     call difstat(nstep,tstep,orb(1,2),orb(1,1),stat(1,1))
      call pgsci(2)
      call difstat(nstep,tstep,orb(1,2),orb(1,3),stat(1,2))
      call pgsci(3)
      call difstat(nstep,tstep,orb(1,2),orb(1,4),stat(1,3))
      call pgsci(4)
      call difstat(nstep,tstep,orb(1,2),orb(1,5),stat(1,4))
      call pgslw(1)
      call pgend

      if (ascout.ne.' ') then
	 nday=(time1-time0)/86400
	 open (10,file=ascout)
	 do i=1,nday
	    call strf1985(string,'%y%m%d',time0+(i-1)*86400)
	    write (10,666) string,(stat(i,j),j=2,4)
	 enddo
      endif
666   format (a,4(1x,f7.2))

      end

      subroutine difstat(n,dt,orb1,orb2,stat)
      integer n,dt
      real*8 orb1(n),orb2(n),stat(*)
      integer nstep,nday,i,j,k
      real*4 rms,rmsold

      nstep=86400/dt
      nday=n/nstep
      k=0
      rms=1e30
      do i=1,nday
         rmsold=rms
	 rms=0
         do j=1,nstep
	    k=k+1
	    rms=rms+(orb1(k)-orb2(k))**2
	 enddo
	 if (rms/nstep.lt.1e5) rms=sqrt(rms/nstep)*100
	 if (rmsold.gt.1e5) then
	    call pgmove(real(i-1),rms)
	 else
	    call pgdraw(real(i-1),rms)
	 endif
	 call pgdraw(real(i),rms)
	 stat(i)=rms
	 if (rms.le.0d0) rms=-999
      enddo
      end
