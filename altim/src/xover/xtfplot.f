      program xtfplot

      integer*4 nmax
      parameter (nmax=50000)
      character filenm*80,spec*4,arg*80,dev*80
      integer*2 i2(4),flags
      integer*4 i4(6),par(10),sat,i,j,m,iargc,lrun,npar,length,nrec
      real*8 u1,u2,pi,rad,rev/0d0/,d
      integer n(nmax),ci(3:6)/2,3,5,6/,mark(3:6)/847,845,846,841/
      real t(nmax,3:6),h(nmax,3:6),tmin/1e35/,tmax/-1e35/,
     .		hmin/1e35/,hmax/-1e35/
      logical tpreset/.false./,hpreset/.false./,color/.false./
      
      data filenm/' '/,dev/' '/

      pi=4*atan(1d0)
      rad=pi/180

      do i=1,iargc()
	 call getarg(i,arg)
	 if (arg(1:2).eq.'t=') then
	    read (arg(3:),*) tmin,tmax
	    tpreset=.true.
	 else if (arg(1:2).eq.'h=') then
	    read (arg(3:),*) hmin,hmax
	    hpreset=.true.
	 else if (arg(1:4).eq.'dev=') then
	    dev=arg(5:)
	 else if (arg.eq.'-color') then
	    color=.true.
	 else if (arg(1:4).eq.'run=') then
	    read (arg(5:),*) lrun
	 else
	    filenm=arg
	 endif
      enddo

      if (filenm.eq.' ') then
	 write (0,1300)
	 goto 9999
      endif

      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=12)
      read (10,rec=1) spec,nrec,npar
      close (10)
      if (spec.eq.'@XTF') then
	 npar=3
      else if (spec.eq.'@XTE' .or. spec.eq.'@XTB') then
      else
	 call fin('Input file is not XTF')
      endif

      length=34+npar*8
      open (10,file=filenm,form='unformatted',status='old',
     .access='direct',recl=length)

      do j=1,nrec
         read (10,rec=j+1) i2,i4,(par(i),i=1,npar*2),flags
	 sat=i2(2)
	 if (sat.eq.3) then
	    rev=6000.0d0
	 else if (sat.eq.4) then
	    rev=6035.0d0
	 else if (sat.eq.5 .or. sat.eq.6) then
	    rev=6745.759d0
	 else
	    write (6,*) i2(2)
	    call fin('unknown sat')
	 endif
	 u1=i4(2)/1d6*rad
	 u2=u1+(i4(6)-i4(5))/rev*2*pi
	 d=par(1)*(u2-u1)
     .		-par(2)*(cos(u2)-cos(u1))
     .		+par(3)*(sin(u2)-sin(u1))
     .		-par(4)*(cos(2*u2)-cos(2*u1))/2
     .		+par(5)*(sin(2*u2)-sin(2*u1))/2
	 d=d/(u2-u1)/1d6
	 if (flags.lt.256) goto 100
	 m=n(sat)+1
	 n(sat)=m
	 t(m,sat)=i4(5)/86400d0/365.25d0+85
	 h(m,sat)=d
	 if (.not.tpreset) then
	    tmin=min(tmin,t(m,sat))
	    tmax=max(tmax,t(m,sat))
	 endif
	 if (.not.hpreset) then
	    hmin=min(hmin,h(m,sat))
	    hmax=max(hmax,h(m,sat))
	 endif
100      continue
      enddo

      call pgbeg(0,dev,1,1)
      call pgswin(tmin,tmax,hmin,hmax)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      call pglab(' ','mean correction (m)',' ')
      if (color) then
	 call pgsci(ci(3))
         call pgmtxt('B',2.5,0.0,0.0,'\\(791) GEOSAT')
	 call pgsci(ci(4))
         call pgmtxt('B',2.5,0.2,0.0,'\\(791) ERS-1')
	 call pgsci(ci(5))
         call pgmtxt('B',2.5,0.4,0.0,'\\(791) TOPEX')
	 call pgsci(ci(6))
         call pgmtxt('B',2.5,0.6,0.0,'\\(791) POSEIDON')
	 call pgsci(1)
      else
         call pgmtxt('B',2.5,0.0,0.0,'\\(847) GEOSAT')
         call pgmtxt('B',2.5,0.2,0.0,'\\(845) ERS-1')
         call pgmtxt('B',2.5,0.4,0.0,'\\(846) TOPEX')
         call pgmtxt('B',2.5,0.6,0.0,'\\(841) POSEIDON')
      endif
      call pgmtxt('B',2.5,1.0,1.0,'time (years)')
      do sat=3,6
	 if (color) then
	    call pgsci(ci(sat))
	    call pgpt(n(sat),t(1,sat),h(1,sat),-1)
	 else
	    call pgpt(n(sat),t(1,sat),h(1,sat),mark(sat))
	 endif
	 if (lrun.ge.2) then
	    call pgslw(3)
	    call runaver(n(sat),t(1,sat),h(1,sat),lrun)
	    call pgline(n(sat),t(1,sat),h(1,sat))
	    call pgslw(1)
	 endif
      enddo
      call pgsci(1)
      call pgend

1300  format ('XTFPLOT plots the mean orbit corrections per track ',
     .'archived in an XTF file'//
     .'usage: xtfplot [options] xtf-filename'//
     .'where [options] are:'/
     .'dev=dev    use plotdevice dev'/
     .'t=y0,y1    only plot time range y0,y1 (years). ',
     .'Default is entire range of XTF'/
     .'h=h0,h1    use vertical range h0,h1 (meters). ',
     .'Default is entire range of XTF'/
     .'run=lrun   plot running average with length of lrun (tracks)'/
     .'-color     use color mode')
9999  end

      subroutine runaver(n,x,y,l)
      implicit none
      integer n,k,m,l,lrun,lrun1,lrun2
      real x(n),y(n)
      real*8 sum,y0

      lrun=min(l,n)
      lrun1=lrun-1
      lrun2=lrun/2

      sum=y(1)
      do k=2,lrun
	 sum=sum+y(k)
      enddo
      y0=y(1)
      y(1)=sum/lrun

      do k=2,n-lrun1
	 sum=sum-y0+y(k+lrun1)
	 y0=y(k)
	 y(k)=sum/lrun
      enddo

      do k=n,1,-1
	 m=max(1,min(k-lrun2,n-lrun1))
	 y(k)=y(m)
      enddo
      end
