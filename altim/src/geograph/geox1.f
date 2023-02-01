      program geox1

* Generate normal equations from mean crossover height differences.
* Input file type is XGF
*
* Software history
* 10-Apr-1998 -- (9804.0) Created out of geox.f
*  9-Aug-2001 -- (0108.0) Adaption to new geosubs.a
*-----------------------------------------------------------------------
      implicit none
      include "geograph.inc"
      integer*4 i,k,l,ip,jp,n,nsel
      integer*4 iarg,iargc,lnblnk
      integer*4 irec,nrec,xgf4(4)
      integer*2 xgf2(9)
      equivalence (xgf2,xgf4)
      parameter (nsel=3715)
      real*4  scale
      real*8  xovarmin/3d-2/
      real*8  x0,x1,y0,y1,dx,dy
      real*8  a(nsel),b(nsel),atb(nsel),d,w,
     |	x,y,cycles,period,h0
      real*8  ata(nsel*(nsel+1)/2)
      character*8  version/'0108.0'/
      character*80 arg,out/' '/,in/' '/,coefs/' '/
      character*4 spec
      namelist /geograph_nml/ incl,period,cycles,deadband,lmax,test,
     | x0,x1,y0,y1,dx,dy

      arg='/user/altim'
      call checkenv('ALTIM',arg,l)
      arg(l+1:)='/nml/geograph.nml'
      open (7,file=arg,status='old')
      read (7,geograph_nml)
      close (7)
      open (7,file="geograph.nml",status="old",err=2)
      read (7,geograph_nml)
      close (7)
2     continue

* Initialise

      lmax=ndeg

* Start printout

      write (*,600) version
600   format('This is geox1 - version ',a)

* Scan arguments

      do iarg=1,iargc()
	 call getarg(iarg,arg)
         if (arg(1:4).eq.'nml=') then
	    l=lnblnk(arg)
	    write (*,550) 'Namelist       : '//arg(5:l)
	    open (7,file=arg(5:),status='old')
	    read (7,geograph_nml)
	    close (7)
         else if (arg(1:5).eq.'lmax=') then
	    read (arg(6:),*) lmax
	 else if (coefs.eq.' ') then
	    coefs=arg
	 else if (in.eq.' ') then
	    in=arg
	 else if (out.eq.' ') then
	    out=arg
	 endif   
      enddo

* Print usage when not all arguments are given

      if (out.eq.' ') then
	 write (*,1300)
	 goto 9999
      endif

* Initialize with apriori coefficients

      call gravrd(0.0,coefs)

* Convert satellite orbit information

      incl=incl*rad
      n0=(2*pi)/(period*86400d0/cycles)
      wmdot=n0
      ogdot=-nint(period)*n0/cycles
      a0=(gm/n0**2)**(1d0/3d0)
      h0=a0-ae

* Initialize Flmp and Dlmp

      call d_lmp

      i=lnblnk(coefs)
      k=lnblnk(in)
      l=lnblnk(out)
      write (*,610) coefs(:i),in(:k),out(:l),lmax,ipmax
610   format(
     |'Coefficients   : ',a/
     |'Input XGF      : ',a/
     |'NormalEquations: ',a/
     |'Max deg/order  : ',i9/
     |'Nr of coeffs   : ',i9)

* Open input XGF

      open (10,file=in,status='old',form='unformatted',
     |		access='direct',recl=18)
      read (10,rec=1) spec,nrec
      read (10,rec=1) xgf2
      
* If construct normal equations. Remember to first set all coefficients
* to 1.

      call statbar (0,nrec,'Building matrix')
      do i=1,ipmax*(ipmax+1)/2
         ata(i)=0d0
      enddo
      do i=1,ipmax
         atb(i)=0d0
	 cs(i)=1d0
      enddo

* Loop: compute arguments for each valid point.
* The array b contains the geographically anti-correlated orbit
* error per coefficient.

      do irec=1,nrec
	 read (10,rec=irec+1) xgf2
	 y=xgf4(2)/1d6*rad
	 x=xgf4(3)/1d6*rad
	 call geocen(y,h0,y,d)	! convert to geocentric
	 w=max(xovarmin,xgf2(9)/1d3)**2/2
	 w=w/xgf2(2)
	 d=xgf4(4)/1d6/2		! estimate half the xover difference
         call geogrcmp(y,x,a,b)
	 do ip=1,ipmax
	    atb(ip)=atb(ip)+b(ip)*d/w
	    do jp=1,ip
	       n=g(ip)+jp
	       ata(n)=ata(n)+b(ip)*b(jp)/w
	    enddo
	 enddo
	 call statbar(1,irec,' ')
      enddo

* Write matrices

      open (20,file=out,form='unformatted')
      write (20) ipmax
      scale=1
      do i=1,ipmax
         write (20) ics(i),ideg(i),iord(i),solve(i),sig(i),scale
      enddo
      do i=1,ipmax
         write (20) (sngl(ata(g(i)+k)),k=1,i)
      enddo
      write (20) (sngl(atb(i)),i=1,ipmax)
      close (20)

550   format(a)
1300  format(
     |'geox1 - compute normal equations from xover height',
     |' differences'//
     |'usage: geox1 [options] coefs xgf equations'//
     |'with'/
     |' coefs     Apriori normal matrix (only used to extract',
     |' coefficients)'/
     |' xgf       XGF input file name'/
     |' equations Normal equations output file name'//
     |'and  [options ]'/
     |' nml=name  Namelist file name (in addition to geograph.nml)',
     |' lmax=lmax Maximum degree/order of coefficients')
9999  end
