      program geospec
      implicit none

      include "geograph.inc"
      integer l,m,i,j,k,p,ip,numb(-ndeg:ndeg,0:ndeg),ios,nstep/0/
      real*8 amp(2,-ndeg:ndeg,0:ndeg)
      real*8 c,s,ampl,x0,x1,y0,y1,dx,dy,period
      integer cycles,iargc,narg,pgbeg,iarg,nfreq,mfreq
      parameter (mfreq=100000)
      real*8 psi0(-ndeg:ndeg,0:ndeg)
      real x,y,fmax/3.0/,amax/0.0/,step/0.0/,ff(mfreq),hh(mfreq)
      character*80 filenm1/' '/,filenm2/' '/,arg,dev/' '/,title/'-'/,
     | table/' '/
      namelist /geograph_nml/ incl,period,cycles,deadband,lmax,test,
     | x0,x1,y0,y1,dx,dy
      namelist /geospec_nml/ fmax,amax,dev,deadband

      arg='/user/altim'
      call checkenv('ALTIM',arg,l)
      arg(l+1:)='/nml/geograph.nml'
      open (7,file=arg,status='old')
      read (7,geograph_nml)
      close (7)
      arg(l+1:)='/nml/geospec.nml'
      open (7,file=arg,status='old')
      read (7,geospec_nml)
      close (7)
      open (7,file="geograph.nml",status='old',err=2)
      read (7,geograph_nml)
      close (7)
2     continue
      open (7,file="geospec.nml",status='old',err=3)
      read (7,geospec_nml)
      close (7)
3     continue

      narg=iargc()
      do iarg=1,narg
	 call getarg(iarg,arg)
	 if (arg(:2).eq.'-t') then
	    test=.true.
	 else if (arg(:4).eq.'dev=') then
	    dev=arg(5:)
	 else if (arg(:5).eq.'amax=') then
	    read (arg(6:),*,iostat=ios) amax,step,nstep
	 else if (arg(:4).eq.'tab=') then
	    table=arg(5:)
	 else if (arg(:6).eq.'title=') then
	    title=arg(7:)
	 else if (filenm1.eq.' ') then
	    filenm1=arg
	 else
	    filenm2=arg
	 endif
      enddo

      if (filenm1.eq.' ') then
	 write (0,600)
600      format('geospec - plot spectrum of gravity-induced radial',
     |' orbit error'//
     |'syntax: geospec [ options ] field1 field2'//
     |'with [options]:'/
     |'tab=name   : specify file name of tabular output'/
     |'amax=amax  : specify maximum amplitude in plot (cm)'/
     |'dev=dev    : specify plot device'/
     |'title=text : specify plot title')
	 goto 9999
      else if (filenm1.eq.'sigma') then
         call gravrd(1.0,filenm2)
      else if (filenm2.eq.' ') then
	 open (10,file=filenm1)
	 rewind (10)
   10    read (10,610,end=20) m,k,numb(k,m),psi0(k,m),
     |		amp(1,k,m),amp(2,k,m)
	 goto 10
   20	 goto 100
      else
         call gravrd(1.0,filenm1)
         call gravrd(-1.0,filenm2)
      endif

*     deadband=0d0

      incl=incl*rad
      n0=(2*pi)/(period*86400d0/cycles)
      wmdot=n0
      ogdot=-nint(period)*n0/cycles
      a0=(gm/n0**2)**(1d0/3d0)

      call d_lmp

      do ip=1,ipmax
         l=ideg(ip)
         m=iord(ip)
         j=ics(ip)
         if (m/2*2.eq.m) then
            do p=0,l
               k=l-2*p
               i=pnt(l)+m*(l+1)+p
	       psi0(k,m)=((l-2*p)*wmdot+m*ogdot)/n0
               amp(j,k,m)=amp(j,k,m)+cs(ip)*dlmp(i)
               numb(k,m)=numb(k,m)+1
            enddo
         else
            do p=0,l
               k=l-2*p
               i=pnt(l)+m*(l+1)+p
	       psi0(k,m)=((l-2*p)*wmdot+m*ogdot)/n0
               amp(3-j,k,m)=amp(3-j,k,m)-(-1)**j*cs(ip)*dlmp(i)
               numb(k,m)=numb(k,m)+1
            enddo
         endif
      enddo

* Combine frequencies with m=0

      do k=0,lmax
	 amp(1,k,0)=amp(1,k,0)+amp(1,-k,0)
	 amp(2,k,0)=amp(2,k,0)-amp(2,-k,0)
	 numb(k,0)=numb(k,0)+numb(-k,0)
	 numb(-k,0)=0
      enddo

* Invert negative frequencies and write out

      do m=0,lmax
         do k=-lmax,lmax
	    if (psi0(k,m).lt.0) then
	       amp(2,k,m)=-amp(2,k,m)
	       psi0(k,m)=-psi0(k,m)
	    endif
	    if (psi0(k,m).lt.fmax .and. numb(k,m).gt.0) then
	    amp(1,k,m)=amp(1,k,m)*100
	    amp(2,k,m)=amp(2,k,m)*100
	    ampl=sqrt(amp(1,k,m)**2+amp(2,k,m)**2)
	    if (abs(psi0(k,m)-1).gt.deadband .and.
     |          psi0(k,m).gt.deadband)
     |         write (*,610) m,k,numb(k,m),
     |         psi0(k,m),amp(1,k,m),amp(2,k,m),ampl
	    endif
	 enddo
      enddo

* Determine spectral lines

  100 continue
      nfreq=0
      do m=0,lmax
         do k=-lmax,lmax
	    if (psi0(k,m).lt.fmax .and. numb(k,m).gt.0
     |		.and. abs(psi0(k,m)-1).gt.deadband .and.
     |		psi0(k,m).gt.deadband) then
	       c=amp(1,k,m)
	       s=amp(2,k,m)
	       x=psi0(k,m)
	       y=sqrt(c*c+s*s)
	       if (y.gt.0 .and. y.lt.1e30) then
	          nfreq=nfreq+1
		  if (nfreq.gt.mfreq) call fin('too many frequencies')
	          ff(nfreq)=x
	          hh(nfreq)=y
	       endif
	    endif
	 enddo
      enddo

* Write out table if requested

      if (table.ne.' ') then
	 open (20,file=table)
	 write (20,620)
         do m=1,nfreq
	    write (20,'(2f8.4)') ff(m),hh(m)
	 enddo
      endif
      if (dev.eq.' ') goto 9999

* Make plot if requested

      if (pgbeg(0,dev,1,1).ne.1)
     |	call fin("can not open plot device")
      m=index(filenm1,' ')-1
      if (title.eq.'-') title=filenm1(:m)//' - '//filenm2
      call pgswin(0.0,fmax,0.0,amax)
      call pgbox('BCNSTI',0.0,0,'BCNSTI',step,nstep)
      call pglab('frequency (cycl/rev)','amplitude (cm)',
     |title)
      do m=1,nfreq
         call pgmove(ff(m),0.0)
	 call pgdraw(ff(m),hh(m))
      enddo
      call pgend

* Formats

610   format (i3,2i4,f10.6,3f8.3)
620   format('# Orbit spectrum'/
     |'# Column 1: frequency (cycles/rev)'/
     |'# Column 2: amplitude (cm)')
9999  end
