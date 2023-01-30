**DISWGRID -- Gaussian distance weighted gridding

      program diswgrid

* Gridding program using Gaussian distance weighting, incl. 
* reinterpolation.
*
* The weight function,
* equals exp(-dist**2/sigma**2) where sigma is constant.
*
* 9510.0 - 22-Oct-1995 - Brand-new version (Remko Scharroo)
* 9510.1 -  1-Nov-1995 - Include minimum weight
* 9512.0 -  7-Dec-1995 - Bug removed that prevented writing of reinterpolated
*			last buffer of the XGF. bad= flag introduced.
* 9711.0 - 13-Nov-1997 - Added temporal weighting. t= flag
* 9802.0 - 24-Feb-1998 - Multi-step gridding
* 9911.0 -  3-Nov-1999 - Added DATEARG function
* 0112.0 -  3-Dec-2001 - Added std= flag
* 0306.0 -  3-Jun-2003 - Proper initialisation of a few variables
* 0308.0 -  6-Aug-2003 - List temporal horizon in output
*----------------------------------------------------------------------
      implicit none

      character*6  version
      parameter    (version='0308.0')
      include 'diswgrid.inc'

      real*8       sum(maxgrd),wgt(maxgrd),gr1int
      integer*4    gridwr8
      integer*4    itg,ixg,iyg,izg
      integer*2    isg
      character*18 xgf
      equivalence  (xgf(1:4),itg),(xgf(5:8),iyg),(xgf(9:12),ixg),
     |             (xgf(13:16),izg),(xgf(17:18),isg)
      integer*4    nused,wrap/0/
      integer*4    i,j,l1,l2,k
      integer*4    kx,kxlo,kxhi,ky,kylo,kyhi
      character    fspec*4,xgfin*80/' '/,gridut*80/' '/,text*13,arg*80,
     |		   addgridnm*80/'-'/,refgridnm*80/' '/
      logical      spherical/.true./,reint/.false./,ref/.false./,
     |             datearg
      real*8	   sigx/0d0/,sigy/0d0/,sigt/-7d0/,hor/2.5d0/,hor2,bad/1d40/
      real*8	   bt/1d40/,wt/0d0/,decr/5d-1/,epoch,dum
      real*8	   bkx,xg,sigkx,x,rx,u/0d0/,ut/0d0/
      real*8	   bky,yg,sigky,y,ry,v/0d0/,vt/0d0/
      real*8	   z,zmin,zmax,zmean,zrms,r,w,s,minwgt/1d-2/
      real*8       smin/1d-3/,smax/1d30/
      integer	   invalid,nrec,iargc,iter/1/,niter/1/

      integer	   openf,readf,ios
      real*8	   pi,rad

* Initialize

      pi=4*atan(1d0)
      rad=pi/180
      x0=-180
      x1=180
      dx=0.1
      y0=-90
      y1=90
      dy=0.1

* Read parameters from command line

      do i=1,iargc()
         call getarg(i,arg)
	 if (arg(:4).eq.'lon=') then
	    read(arg(5:),*,iostat=ios) x0,x1,dx
	 else if (arg(:4).eq.'lat=') then
	    read(arg(5:),*,iostat=ios) y0,y1,dy
	 else if (arg(:4).eq.'wgt=') then
	    read(arg(5:),*,iostat=ios) minwgt
	 else if (arg(:4).eq.'hor=') then
	    read(arg(5:),*,iostat=ios) hor
         else if (arg(:4).eq.'bad=') then
	    read(arg(5:),*,iostat=ios) bad
         else if (arg(:4).eq.'std=') then
	    read(arg(5:),*,iostat=ios) smin,smax
	 else if (datearg(arg,epoch,dum,dum)) then
	    l1=index(arg,'=')+1
	    read(arg(l1:),*,iostat=ios) dum,sigt,u,v
	    sigt=abs(sigt*86400)
         else if (arg(:4).eq.'sig=') then
	    read(arg(5:),*,iostat=ios) sigx,sigy
	    if (sigy.eq.0) sigy=sigx
	 else if (arg(:4).eq.'ref=') then
	    refgridnm=arg(5:)
	    ref=.true.
	 else if (arg(:4).eq.'add=') then
	    addgridnm=arg(5:)
	 else if (arg(:5).eq.'iter=') then
	    read(arg(6:),*,iostat=ios) niter,decr
         else if (arg(:2).eq.'-i') then
	    reint=.true.
	 else if (arg(:2).eq.'-r') then
	    spherical=.false.
	 else if (arg(:2).eq.'-w') then
	    wrap=1
         else if (xgfin.eq.' ') then
	    xgfin=arg
         else
	    gridut=arg
         endif
      enddo

      if (xgfin.eq.' ' .or. gridut.eq.' ') goto 1360

      l1=index(xgfin,' ')-1
      l2=index(gridut,' ')-1
      write (*,1000) version,xgfin(:l1),gridut(:l2)
      l1=index(refgridnm,' ')-1
      l2=index(addgridnm,' ')-1
      if (ref) write (*,1005) ' Ref',refgridnm(:l1)
      if (addgridnm.ne.'-') write (*,1005) ' Add',addgridnm(:l2)

* Perform some preparatory computations and tests.

      nx=nint((x1-x0)/dx)+1
      dx=(x1-x0)/(nx-1)
      ny=nint((y1-y0)/dy)+1
      dy=(y1-y0)/(ny-1)

      text='(rectangular)'
      if (spherical) text='(spherical)'

      if (sigx.eq.0) sigx=2*dx
      if (sigy.eq.0) sigy=2*dy

      if (nx.lt.1 .or. ny.lt.1 .or. nx*ny.gt.maxgrd) goto 1302

* Load or set reference grid

      do k=1,nx*ny
         h(k)=0
      enddo
      if (ref) call addgrid(refgridnm,wgt,maxgrd)

100   continue

      write (*,1010) iter,x0,x1,dx,sigx,hor*sigx,text,
     |	y0,y1,dy,sigy,hor*sigy,minwgt,wrap
      call strf1985(xgf,'%y%m%d %H:%M:%S',nint(epoch))
      if (sigt.gt.0) write (*,1020)
     |  xgf,sigt/86400,sigt*hor/86400,u,v

* Transform the units of sigmas and boundaries
* from degrees to grid distances

      sigkx=sigx/dx
      sigky=sigy/dy
      bkx=hor*sigkx
      bky=hor*sigky
      hor2=hor*hor

* Prepare variables for temporal selection and weighting

      if (sigt.gt.0) bt=hor*sigt

* Open input XGF file.

      fd=openf(xgfin,'r')
      if (fd.lt.0) goto 1307
      ios=readf(fd,4,fspec)
      rec1=0
      if (fspec.ne.'@XGF') goto 1350
      ios=readf(fd,4,nrec)

* Prepare arrays for grid; set counters

      do k=1,nx*ny
         sum(k)=0
	 wgt(k)=0
      enddo
      nused=0
      zmin=+1d40
      zmax=-1d40
      zmean=0
      zrms=0
      call statbar(0,nrec,'Gridding')

* Read and process XGF file

      do i=1,nrec
         call xgfrd(i,xgf)
	 if (isg.lt.0) goto 220

* Store lat, lon, sigma in units of degrees and metres

         x=ixg*1d-6
         y=iyg*1d-6
	 s=isg*1d-3

* Check time selection and compute temporal weight and horizontal
* shift for this data point

         if (sigt.le.0) then
	 else if (abs(itg-epoch).gt.bt) then
	    goto 220
	 else
	    ut=(itg-epoch)*u/111d5
	    vt=(itg-epoch)*v/111d5/cos(y*rad)
	    wt=((itg-epoch)/sigt)**2
	 endif

* Limit standard deviation to given limits

	 s=max(smin,min(s,smax))

* Transform x-sigma in case of spherical distance weighting

         if (spherical) then
            sigkx=sigx/dx/cos(y*rad)
            bkx=hor*sigkx
	 endif

* Compute test indices for this data point;
* apply the advection terms here.

         yg=(y-vt-y0)/dy+1
         kylo=max(nint(yg-bky+0.5), 1)
         kyhi=min(nint(yg+bky-0.5),ny)
	 if (kyhi.lt.kylo) goto 220

         do j=-wrap,wrap
	 xg=(x-ut-x0+j*360d0)/dx+1
	 kxlo=max(nint(xg-bkx+0.5), 1)
	 kxhi=min(nint(xg+bkx-0.5),nx)
	 if (kxhi.lt.kxlo) goto 210

* Interpolate from reference if required

         if (ref) then
	    r=gr1int(h,xg,yg,nx,ny,nx)
	    if (r.gt.1d20) goto 210
	    z=izg*1d-6-r
	 else
	    z=izg*1d-6
	 endif

* Keep residual statistics

         nused=nused+1
         zmax=max(zmax,z)
	 zmin=min(zmin,z)
         zmean=zmean+z
         zrms=zrms+z*z

* Spread data point over neighbouring grid points

         do ky=kylo,kyhi
            ry=(yg-ky)/sigky
            do kx=kxlo,kxhi
               rx=(xg-kx)/sigkx
               r=rx*rx+ry*ry
               if (r.le.hor2) then
                  w=exp(-r-wt)/s
                  k=kx+(ky-1)*nx
                  wgt(k)=wgt(k)+w
                  sum(k)=sum(k)+w*z
               endif
	    enddo
         enddo

210	 continue
	 enddo

220	 call statbar(1,i,' ')
      enddo

      zmean=zmean/nused
      zrms=sqrt(zrms/nused-zmean**2)
      write (*,1070) nrec,zmin,zmean,nused,zmax,zrms

* Postprocess grid

      if (nused.eq.0) write (*,550) 'no data to process.'
      invalid=0
      nused=0
      zmin=+1d40
      zmax=-1d40
      zmean=0
      zrms=0

      do k=1,nx*ny
         w=wgt(k)
         if (w.gt.minwgt) then
	    z=sum(k)/w
	    zmin=min(zmin,z)
	    zmax=max(zmax,z)
	    zmean=zmean+z
	    zrms=zrms+z*z
	    nused=nused+1
         else
            z=bad
            invalid=invalid+1
         endif
	 h(k)=h(k)+z
      enddo

      if (ref) then

      zmean=zmean/nused
      zrms=sqrt(zrms/nused-zmean**2)
      write (*,1075) 'Incremental grid statistics:',
     |		nused,zmin,zmean,invalid,zmax,zrms

      invalid=0
      nused=0
      zmin=+1d40
      zmax=-1d40
      zmean=0
      zrms=0

      do k=1,nx*ny
         z=h(k)
         if (z.lt.1d20) then
	    zmin=min(zmin,z)
	    zmax=max(zmax,z)
	    zmean=zmean+z
	    zrms=zrms+z*z
	    nused=nused+1
         else
            invalid=invalid+1
         endif
      enddo

      endif

      zmean=zmean/nused
      zrms=sqrt(zrms/nused-zmean**2)
      write (*,1075) 'Total grid statistics:',
     |		nused,zmin,zmean,invalid,zmax,zrms

* Reloop in case of more iterations

      if (iter.lt.niter) then
         ref=.true.
	 bad=0d0
	 sigx=sigx*decr
	 sigy=sigy*decr
	 iter=iter+1
	 goto 100
      endif

* Add grid to final result if required

      if (addgridnm.ne.'-') call addgrid(addgridnm,wgt,maxgrd)

* Now start output of grid

      if (gridwr8(gridut,nx,ny,h,nx,x0,x1,y0,y1).ne.0)
     |		write (*,550) 'grid probably unusable'
      
      if (.not.reint) goto 9999

* Initialization for re-interpolation

      zmin=+1d40
      zmax=-1d40
      zmean=0
      zrms=0
      rec1=0
      nused=0
      call statbar (0,nrec,'Reinterpolating')

* Reinterpolation pass

*     call seekf(fd,18,0)
      do i=1,nrec
	 call xgfrd(i,xgf)

         x=ixg*1d-6
         y=iyg*1d-6
         z=izg*1d-6

         xg=(x-x0)/dx+1
         yg=(y-y0)/dy+1
         r=gr1int(sum,xg,yg,nx,ny,nx)
	 if (r.lt.1e20) then
	    z=z-r
            izg=nint(z*1d6)
	 else if (isg.gt.0) then
	    isg=-isg
	 endif

	 call xgfwr(i,xgf)

	 if (isg.gt.0) then
	    nused=nused+1
            zmax=max(zmax,z)
            zmin=min(zmin,z)
            zmean=zmean+z
            zrms=zrms+z**2
	 endif

	 call statbar(1,i,' ')
      enddo

      zmean=zmean/nused
      zrms=sqrt(zrms/nused-zmean**2)
      write (*,1070) nrec,zmin,zmean,nused,zmax,zrms
      goto 9999

* Formats

  550 format (a)
 1000 format ('This is DISWGRID, version ',a//
     |' XGF : ',a/'Grid : ',a)
 1005 format (a,' : ',a)
 1010 format (/32('*'),' Iteration',i3,1x,32('*')//
     |11x,'  Minimum  Maximum Cellsize    Sigma  Horizon'/
     |'Longitude :',5f9.3,2x,a/'Latitude  :',5f9.3/
     |'Min weight:',f9.3,'/meter'/
     |'Wrapping  :',i9,' (0=off, 1=on)')
 1020 format (
     |'Epoch     :',3x,a15,9x,2f9.3/
     |'Speed u,v :',2f9.3,' cm/s')
 1070 format (/'Measurement residuals:'/
     |i8,' data points read       min:',f9.3,5x,'mean:',f9.3/
     |i8,' data points used       max:',f9.3,4x,'sigma:',f9.3)
 1075 format(/a/
     |i8,' valid grid points      min:',f9.3,5x,'mean:',f9.3/
     |i8,' invalid grid points    max:',f9.3,4x,'sigma:',f9.3)

* Errors

 1302 write (*,550) 'illegal grid dimension.'
      goto 9999

 1307 write (*,550) 'error opening XGF-file.'
      goto 9999

 1350 write (*,550) 'input data not in XGF format.'
      goto 9999

 1360 write (*,1361) version
 1361 format('DISWGRID, version ',a,': Spatio-temporal gridding'//
     |'usage: diswgrid [ options ] xgf-filename grid-filename'//
     |'where [ options ] are:'/
     |'lon=x0,x1,dx   : longitude boudaries and grid cell size (deg)'/
     |'lat=y0,y1,dy   : latitude  boudaries and grid cell size (deg)'/
     |'sig=sigx[,sigy]: sigma in Gaussian weighting function (deg)',
     |' (def: 2*dx,2*dy)'/
     |'t=t,sigt,u,v   : weight data to epoch with time sigma (days)',
     |' and advection (cm/s)'/
     |17x,'(default is no temporal weighting nor advection)'/
     |'             ... or use mjd=, doy=, ymd=, sec='/
     |'hor=hor        : horizon (-) (def: 2.5)'/
     |'bad=bad        : specify bad value in 1st iter. (def: 1d30)'/
     |'                 (for later iterations bad=0)'/
     |'wgt=wgt        : specify minimum weight for determined points',
     |' (def: 0.01/meter)'/
     |'std=smin,smax  : specify limits on measurement stddev (m)'/
     |'ref=gridnm     : start with grid as first guess | need not be'/
     |'add=gridnm     : add grid after final result    | same size'/
     |'iter=n[,decr]  : iterate with decreasing sigma (def: 1,0.5)'/
     |'-wrap          : wrap data along date boundary (def: off)'/
     |'-rect          : rectangular distance (def: spherical)'/
     |'-int           : re-interpolate (def: off)')

 9999 call closef(fd)
      end
*********************************************************************
      subroutine xgfrd(rec,xgf)
      integer rec,ios,seekf,readf
      character*18 xgf
      include "diswgrid.inc"

      if (rec.gt.rec1) then
         rec0=rec-1
         ios=seekf(fd,(rec0+1)*18,0)
	 ios=readf(fd,buflen*18,xgfbuf)
	 rec1=rec0+ios/18
      endif
      xgf=xgfbuf(rec-rec0)
      end

      subroutine xgfwr(rec,xgf)
      integer rec,seekf,writef,ios
      character*18 xgf
      include "diswgrid.inc"

      xgfbuf(rec-rec0)=xgf
      if (rec.eq.rec1) then
	 ios=seekf(fd,(rec0+1)*18,0)
	 ios=writef(fd,(rec1-rec0)*18,xgfbuf)
      endif
      end
*********************************************************************
      subroutine addgrid(gridnm,grid,mxp)
      character*80 gridnm
      integer mxp,nxp,nyp,ios,kx,ky,k,gridrd8
      real*8 grid(mxp),x0p,x1p,dxp,y0p,y1p,dyp,z0p,z1p,xp,yp,gr1int
      include 'diswgrid.inc'

      nxp=0
      nyp=mxp
      ios=gridrd8(gridnm,nxp,nyp,grid,x0p,x1p,y0p,y1p,z0p,z1p)
      if (ios.ne.0) stop 'addgrid: error loading grid'
      dxp=(x1p-x0p)/(nxp-1)
      dyp=(y1p-y0p)/(nyp-1)
      do ky=1,ny
	 do kx=1,nx
	    k=(ky-1)*nx+kx
	    xp=(x0+(kx-1)*dx-x0p)/dxp+1
	    yp=(y0+(ky-1)*dy-y0p)/dyp+1
	    h(k)=h(k)+gr1int(grid,xp,yp,nxp,nyp,nxp)
	 enddo
      enddo
      end
