      program geox2

* Estimate gravity field coefficients to minimize geographically correlated
* mean crossover height difference.
* Input file type is XGF
*
* Software history
* 11-Apr-1998 -- (9804.0) Created out of geox.f.
* 14-Apr-1998 -- (9804.1) Also to be used to verify gravity field
* differences
*  9-Aug-2001 -- (0108.0) Adaption to new geosubs.a
*-----------------------------------------------------------------------
      implicit none
      include "geograph.inc"
      integer*4 i,k,l,m,ip,nsel
      integer*4 iarg,iargc,lnblnk
      integer*4 irec,nrec,xgf4(4)
      integer*2 xgf2(9)
      equivalence (xgf2,xgf4)
      parameter (nsel=npar)
      real*8  x0,x1,y0,y1,dx,dy
      integer nx,ny,mx,kx,ky
      parameter (mx=361*181)
      real*8  grid_s(mx),grid_c(mx)
      real*8  a(nsel),d,z,r,w,x,y,cycles,period,h0
      real*8  zmin,zmax,zmean,zrms
      real*8  dmin,dmax,dmean,drms
      real*8  rmin,rmax,rmean,rrms
      real*8  wmin,wmax,wmean,wrms
      real*4  ata(nsel*(nsel+1)/2),atb(nsel),xvect(nsel)
      equivalence (ata,covar),(atb,rvect)
      real*4  fact
      integer neq/0/
      character*6  version/'0108.0'/
      character*80 arg,out/' '/,in/' '/,ceq(10)
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
      
* Print usage notes when no arguments are given

      if (iargc().eq.0) then
	 write (*,1300)
	 goto 9999
      endif

* Start printout

      write (*,600) version
600   format('This is geox2 - version ',a/)

* Scan arguments

      do iarg=1,iargc()
	 call getarg(iarg,arg)
         if (arg(1:4).eq.'nml=') then
	    l=lnblnk(arg)
	    write (*,550) 'Namelist       : '//arg(5:l)
	    open (7,file=arg(5:))
	    read (7,geograph_nml)
	    close (7)
	 else if (arg(1:5).eq.'lmax=') then
	    read (arg(6:),*) lmax
	 else if (arg(1:3).eq.'in=') then
	    in=arg(4:)
	    l=lnblnk(arg)
	    write (*,550) 'Input XGF file : '//arg(4:l)
	 else if (arg(1:4).eq.'out=') then
	    out=arg(5:)
	    l=lnblnk(arg)
	    write (*,550) 'Output directory: '//arg(5:l)
         else
	    l=lnblnk(arg)
	    write (*,550) 'Normal equations: '//arg(:l)
	    m=index(arg,'*')
	    if (m.gt.0) then
	       read (arg(1:m-1),*) fact
	    else
	       fact=1
	    endif      
	    call gravrd(fact,arg(m+1:))
	    neq=neq+1
	    ceq(neq)=arg(m+1:)
	 endif
      enddo

* Do some checks

      nx=nint((x1-x0)/dx)+1
      ny=nint((y1-y0)/dy)+1
      if (nx*ny.gt.mx) call fin('too many gridpoints')

* Convert satellite orbit information

      incl=incl*rad
      n0=(2*pi)/(period*86400d0/cycles)
      wmdot=n0
      ogdot=-nint(period)*n0/cycles
      a0=(gm/n0**2)**(1d0/3d0)
      h0=a0-ae

* Initialize Flmp and Dlmp

      call d_lmp

* Solve normal equation

      l=lnblnk(out)
      if (atb(1).ne.0) then
*      do i=1,ipmax*(ipmax+1)/2
*         ata(i)=covar(i)
*      enddo
*      do i=1,ipmax
*         atb(i)=rvect(i)
*      enddo	 	 
*      call lincho(ata,g,atb,a,b,ipmax,d,i)
      call spptrf('U',ipmax,ata,i)
      if (i.ne.0) call fin('error return from SPPTRF')
      call spptrs('U',ipmax,1,ata,atb,ipmax,i)
      if (i.ne.0) call fin('error return from SPPTRS')

* Print out the solution

      open (13,file=out(:l)//'/Coefficients')
      call statnull(dmin,dmax,dmean,drms)
      do ip=1,ipmax
	 d=atb(ip)/sig(ip)
	 xvect(ip)=atb(ip)
	 write (13,700) ip,ideg(ip),iord(ip),ics(ip),atb(ip),d
	 cs(ip)=atb(ip)
	 if (ip.le.ipmax) call statup(d,dmin,dmax,dmean,drms)
      enddo
      dmean=dmean/ipmax
      drms=sqrt(drms/ipmax)
      write (*,710) ipmax,dmin,dmax,dmean,drms
      close (13)
      else
         do i=1,ipmax
	    cs(i)=-cs(i)
	 enddo
      endif

* Open input XGF

      if (in.ne.' ') then
      open (10,file=in,status='old',form='unformatted',
     |		access='direct',recl=18)
      read (10,rec=1) spec,nrec

* Write out the computed anti-correlated orbit error and
* residuals. Keep statistics.

      read (10,rec=1) xgf2
      open (11,file=out(:l)//'/Correction.xgf',form='unformatted',
     |		access='direct',recl=18)
      write (11,rec=1) xgf2
      open (12,file=out(:l)//'/Residuals.xgf',form='unformatted',
     |		access='direct',recl=18)
      write (12,rec=1) xgf2
      call statnull(dmin,dmax,dmean,drms)
      call statnull(zmin,zmax,zmean,zrms)
      call statnull(rmin,rmax,rmean,rrms)
      call statnull(wmin,wmax,wmean,wrms)
      do irec=1,nrec
	 read (10,rec=irec+1) xgf2
	 y=xgf4(2)/1d6*rad
	 x=xgf4(3)/1d6*rad
	 z=xgf4(4)/1d6
	 call geocen(y,h0,y,d)		! convert to geocentric
         call geograph(y,x,a,d)
	 d=2*d
	 r=z-d
	 xgf4(4)=nint(d*1d6)
	 write (11,rec=irec+1) xgf2
	 xgf4(4)=nint(r*1d6)
	 write (12,rec=irec+1) xgf2
	 call statup(z,zmin,zmax,zmean,zrms)
	 call statup(d,dmin,dmax,dmean,drms)
	 call statup(r,rmin,rmax,rmean,rrms)
	 call statup(w,wmin,wmax,wmean,wrms)
      enddo
      close (10)
      close (11)
      close (12)

      zmean=zmean/nrec
      zrms =sqrt(zrms/nrec-zmean**2)
      dmean=dmean/nrec
      drms =sqrt(drms/nrec-dmean**2)
      rmean=rmean/nrec
      rrms =sqrt(rrms/nrec-rmean**2)
      wmean=wmean/nrec
      wrms =sqrt(wrms/nrec-wmean**2)

      write(*,550)'Measurement statistics (cm)'
      write(*,720)'Input XGF',nrec,zmin*100,zmax*100,zmean*100,zrms*100
      write(*,720)'      2*S',nrec,dmin*100,dmax*100,dmean*100,drms*100
      write(*,720)'Residuals',nrec,rmin*100,rmax*100,rmean*100,rrms*100
      endif

* Grid the computed anti-correlated orbit error.

      call statnull(dmin,dmax,dmean,drms)
      call statnull(zmin,zmax,zmean,zrms)
      call statnull(wmin,wmax,wmean,wrms)
      k=0
      do ky=1,ny
         y=(y0+(ky-1)*dy)*rad
	 call geocen(y,h0,y,d)		! convert to geocentric
	 call statup(w,wmin,wmax,wmean,wrms)
         do kx=1,nx
	    k=k+1
            x=(x0+(kx-1)*dx)*rad
            call geograph(y,x,z,d)
	    d=2*d
            grid_c(k)=z
            grid_s(k)=d
	    call statup(d,dmin,dmax,dmean,drms)
	    call statup(z,zmin,zmax,zmean,zrms)
         enddo
      enddo
      call gridwr8(out(:l)//'/Geocor_c.grd',nx,ny,grid_c,nx,x0,x1,y0,y1)
      call gridwr8(
     |out(:l)//'/Geocor_2s.grd',nx,ny,grid_s,nx,x0,x1,y0,y1)
      nrec=nx*ny
      zmean=zmean/nrec
      zrms =sqrt(zrms/nrec-zmean**2)
      dmean=dmean/nrec
      drms =sqrt(drms/nrec-dmean**2)
      wmean=wmean/ny
      wrms =sqrt(wrms/ny-wmean**2)
      write(*,550)'Grid statistics (cm)'
      write(*,720)'        C',nrec,zmin*100,zmax*100,zmean*100,zrms*100
      write(*,720)'      2*S',nrec,dmin*100,dmax*100,dmean*100,drms*100

* Compute reduction of variances ( bTb + xT(ATAx - 2ATb) )

      if (atb(i).ne.0) then
      write (*,"(/'Reduction of variances (bTb)')")
      do i=1,neq
         do k=1,ipmax*(ipmax+1)/2
	    ata(k)=0
	 enddo
         do k=1,ipmax
	    atb(k)=0
	 enddo
         call gravrd(1.0,ceq(i))
	 call sspmv('U',ipmax,-1.0,ata,xvect,1,2.0,atb,1)
	 d=0
	 do k=1,ipmax
	    d=d+xvect(k)*atb(k)
	 enddo
	 l=lnblnk(ceq(i))
	 write (*,730) d,ceq(i)(:l)
      enddo
      endif

550   format (a)
700   format (i4,3i3,d18.9,f14.9)
710   format (/'Gravity adjustment:',i6,5f12.6/)
720   format (a9,' : ',i9,4f12.3)
730   format (d18.9,' : ',a)

1300  format(
     |'geox2 - compute graviy field correction from normal equations'//
     |'usage: geox2 [ options ] fact*normaleqn',
     |' [ fact*normaleqn ... ] out'//
     |'with'/
     |' fact      Factor to multiply normal-equations'/
     |' normaleqn Normal-equations file name'/
     |' out       Output directory'//
     |'and with [options]'/
     |' nml=nml   Use namelist nml in combination with geograph.nml'/
     |' lmax=lmax Maximum degree/order'/
     |' out=dir   Output will be stored in directory dir (def:geox2)')
   
9999  end

      subroutine statup(z,zmin,zmax,zmean,zrms)
      real*8 z,zmin,zmax,zmean,zrms
      zmin =min(zmin,z)
      zmax =max(zmax,z)
      zmean=zmean+z
      zrms =zrms+z*z
      end

      subroutine statnull(zmin,zmax,zmean,zrms)
      real*8 zmin,zmax,zmean,zrms
      zmin =+1d40
      zmax =-1d40
      zmean=0
      zrms =0
      end
