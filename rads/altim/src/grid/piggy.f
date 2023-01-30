**PIGGY -- Spatio-temporal gridding program
*+
      program piggy

* This is a program to do spatio-temporal gridding of along-track data
*-
* 22-Feb-1998 - Remko Scharroo
*  3-Nov-1999 - Added DATEARG function
*-----------------------------------------------------------------------
      implicit none
      include 'piggy.inc'
      real*8    dh(mx),sum(mx),rmax,tmax,dum,
     |          horizon,tsigma,dr2,dt,rweight,tweight,weight,
     |          glat,glon,mint,epgrid/0d0/,
     |          rmsgr,mingr,maxgr,rmsh,minh,maxh,diff,rlon,rlat,
     |          rad,pi,height,gr3int,rstep,rsigma,
     |		rsigma0/1d0/,rsigma1/1d0/,rsigma2,hor2,coslat
      integer*4   iarg,iargc,iobs,nobs,ibad,k,
     |          ios,fd,openf,closef,readf,iread,igood,
     |          ix,ix0,ix1,iy,iy0,iy1,hx,hy,cx,cy,itrack,
     |          epobs,minep,maxep,xglen,gridwr8,
     |          iter,nxp,nyp,gridrd8,niter/1/,nfile/0/
      character*80 gridnm,arg,xgfnm,ftype*4,
     |		refgridnm/'-'/,addgridnm/'-'/,filenm(2)
      logical   colt/.false./,datearg
      integer*2 xgf2(9)
      integer*4 xgf4(4)
      equivalence (xgf2,xgf4)
*
* Initialise some constants
*
      pi=4d0*datan(1d0)
      rad=pi/180
      tsigma=7
      horizon=2.5
      iter=0
*
* Get command line arguments
*
      do iarg=1,iargc()
	 call getarg(iarg,arg)
	 if (datearg(arg,epgrid,dum,dum)) then
	    k=index(arg,'=')+1
	    read (arg(k:),*,iostat=ios) dum,tsigma
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*,iostat=ios) lat0,lat1,dlat
	 else if (arg(:4).eq.'lon=') then
	    read (arg(5:),*,iostat=ios) lon0,lon1,dlon
	 else if (arg(:4).eq.'hor=') then
	    read (arg(5:),*) horizon
	 else if (arg(:7).eq.'rsigma=') then
	    read (arg(8:),*,iostat=ios) rsigma0,rsigma1,niter
	    if (rsigma1.eq.0) rsigma1=rsigma0
	 else if (arg(:4).eq.'-col') then
	    colt=.true.
	 else if (arg(:4).eq.'ref=') then
	    refgridnm=arg(5:)
	 else if (arg(:4).eq.'add=') then
	    addgridnm=arg(5:)
	 else
	    nfile=nfile+1
	    filenm(nfile)=arg
	 endif
      enddo

      if (nfile.ne.2) then
         write (0,1300)
	 goto 9999
      endif
1300  format('PIGGY -- Spatio-temporal gridding program'//
     |'syntax: piggy [options] xgfname(s) gridname'//
     |'with required arguments:'/
     |' xgfname(s)           : Names of the input XGF files'/
     |' gridname             : Name of the output grid file'//
     |'and [options] are :'/
     |' t=epoch[,tsigma]     : Grid epoch and time sigma in days',
     |' (def:850101.0,7)'/
     |'                    ... or use mjd=, doy=, ymd=, sec='/
     |' lat=lat0,lat1[,dlat] : Latitude boudaries and cellsize',
     |' (def:-90,90,1)'/
     |' lon=lon0,lon1[,dlon] : Longitude boudaries and cellsize',
     |' (def:-180,180,1)'/
     |' rsigma=sig0[,sig1,n] : Largest/smallest spatial sigma and'/
     |'                        number of gridding steps (def:1,1,1)'/
     |' hor=horizon          : Spatio-temporal horizon (def:2.5)'/
     |' ref=gridname         : Start with first guess (def:none)'/
     |' add=gridname         : Add grid to final result (def:none)'/
     |' -col                 : Input type is colinear XGF')
*
* Compute actual gridsizes from area boundaries and gridpoints per degree
*
      nx=nint((lon1-lon0)/dlon)+1
      ny=nint((lat1-lat0)/dlat)+1
      dlon=(lon1-lon0)/(nx-1)
      dlat=(lat1-lat0)/(ny-1)
*
* Compute the date and the time frame for the data
* Data outside this time frame is to be ignored
*
      tsigma=tsigma*864d2
      tmax=tsigma*horizon
      minep=nint(epgrid-tmax)
      maxep=nint(epgrid+tmax)
*
* Initialize grid
*
      do k=1,nx*ny
	 h(k)=0
      enddo
      if (refgridnm.ne.'-') call addgrid(refgridnm,sum,mx)
*
* Start loop
*
700   continue
      iter=iter+1
      if (iter.eq.1) then
         rsigma=rsigma0
      else
         rsigma=exp(log(rsigma0)+
     |	(log(rsigma1)-log(rsigma0))/(niter-1)*(iter-1))
      endif
      do k=1,nx*ny
	 dh(k)=0d0
	 sum(k)=0d0
      enddo
      rsigma2=rsigma**2
      hor2=horizon**2
      rmax=rsigma*horizon
*
* Open input file, read and check header, and determine filetype
*
      fd=openf(filenm(1),'r')
      if (fd.lt.0) then
         write (0,550) 'piggy: error opening XGF file '//filenm(1)
550      format (a)
	 stop
      endif
      ios=readf(fd,4,ftype)
      ios=readf(fd,4,nobs)
      ios=readf(fd,10,xgf2)
      xglen = 18
*
* Initialisation for each loop
*
      iread=0
      igood=0
      mint=1d30
      rmsh=0d0
      minh=1d30
      maxh=-1d30
      write (arg,'("Iter ",i2,": Scale",f5.2)') iter,rmax
      if (iter.eq.1) then
	 call statbar(0,nobs,arg)
      else
	 call statbar(-1,nobs,arg)
      endif

* Begin first loop, read all measurements and add them to the
* grids if they are within the specified timeframe

      do iobs=1,nobs
	 ios=readf(fd,xglen,xgf2)
	 if (ios.ne.xglen) stop 'piggy: error reading XGF'

         if (colt) then
	    if (xgf2(9).gt.0) goto 100
         else
	    if (xgf2(9).lt.0) goto 100
         endif

	 epobs=xgf4(1)
	 if (epobs.lt.minep) goto 100
	 if (epobs.gt.maxep) goto 100

	 rlat=xgf4(2)/1d6
	 rlon=xgf4(3)/1d6
	 height=xgf4(4)/1d6

	 igood=igood+1
	 iread=iread+1
*
* Statistics of good difference from mean
*
	 rmsh=rmsh+height**2
	 minh=min(minh,height)
	 maxh=max(maxh,height)
*
* Determine "window" around the altimeter data location
*
	 coslat=dcos(rlat*rad)
	 cx=nint((rlon-lon0)/dlon+1)
	 cy=nint((rlat-lat0)/dlat+1)
	 hy=int(rmax/dlat)+1
	 hx=int(rmax/dlon/coslat)+1
	 ix0=max(1,cx-hx)
	 iy0=max(1,cy-hy)
	 ix1=min(nx,cx+hx)
	 iy1=min(ny,cy+hy)
*
* Determine weight function in temporal sense
*
	 tweight=exp(-((epobs-epgrid)/tsigma)**2)
	 mint=min(mint,dt)
*
* Add difference from current result to surrounding grid points
*
	 glon=(rlon-lon0)/dlon+1
	 glat=(rlat-lat0)/dlat+1
	 height=height-gr3int(h,glon,glat,nx,ny,nx)

	 do iy=iy0,iy1
	    glat=lat0+(iy-1)*dlat
	    do ix=ix0,ix1
	       glon=lon0+(ix-1)*dlon
	       dr2=(((glon-rlon)*coslat)**2+(glat-rlat)**2)/rsigma2
	       if (dr2.le.hor2) then
		  rweight=exp(-dr2)
		  weight=rweight*tweight
	          k=(iy-1)*nx+ix
	          dh (k)=dh (k)+height*weight
	          sum(k)=sum(k)+weight
	       endif
	    enddo
	 enddo
*
* Go to next altimeter observation
*

100      call statbar(1,iobs,' ')
      enddo

      ios=closef(fd)
*
* Compute and print the data statistics
*
      rmsh=sqrt(rmsh/igood)
      if (iter.eq.niter) then
      write (6,2030) iread,iread-igood,mint/864d2,
     |               rmsh*100,minh*100,maxh*100
      endif
*
* Fill grid
*
      ibad=0
      igood=0
      rmsgr=0d0
      mingr=1d30
      maxgr=-1d30
      do iy=1,ny
         do ix=1,nx
	    k=(iy-1)*nx+ix
	    if (sum(k).gt.0d0) then
	       dh(k)=dh(k)/sum(k)
	       igood=igood+1
	       h(k)=h(k)+dh(k)
	       mingr=min(mingr,h(k))
	       maxgr=max(maxgr,h(k))
	       rmsgr=rmsgr+h(k)**2
	    else
	       dh(k)=1d30
	       ibad=ibad+1
	    endif
	 enddo
      enddo
*
* Compute and print the data statistics
*
      rmsgr=sqrt(rmsgr/igood)
      if (iter.eq.niter) then
      write (6,2040) igood,ibad,
     |               rmsgr*100,mingr*100,maxgr*100
      endif
      if (iter.lt.niter) goto 700
*
* Add final grid if required
*
      if (addgridnm.ne.'-') call addgrid(addgridnm,sum,mx)
*
* Write resulting grid
*
  710 ios=gridwr8(filenm(2),nx,ny,h,nx,lon0,lon1,lat0,lat1)
      if (ios.ne.0) go to 1330
      goto 9999

 800  write (6,*) 'Error opening plusgrid'
      goto 9999

 1330 write (6,*)'Error writing grid'
      goto 9999

 2030 format (/'Data statistics'/
     |        '---------------'/
     |        'Observation read  :',i8/
     |        'Bad observations  :',i8/
     |        'Minimum time gap  :',f8.1,' days'/
     |        'RMS               :',f8.1,' cm'/
     |        'Minimum           :',f8.1,' cm'/
     |        'Maximum           :',f8.1,' cm')
 2040 format (/'Grid statistics'/
     |        '---------------'/
     |        'Good gridpoints   :',i8/
     |        'Bad gridpoints    :',i8/
     |        'RMS               :',f8.1,' cm'/
     |        'Minimum           :',f8.1,' cm'/
     |        'Maximum           :',f8.1,' cm')
 3000 format ('Iteration number',(1x,i4,1x),
     |           'Spatial scale',(1x,f4.2,1x))

 9999 end
*********************************************************************
      subroutine addgrid(gridnm,grid,mxp)
      character*80 gridnm
      integer mxp,nxp,nyp,ios,ix,iy,k,gridrd8
      real*8 grid(mxp),x0p,x1p,dxp,y0p,y1p,dyp,z0p,z1p,xp,yp,gr3int
      include 'piggy.inc'

      nxp=0
      nyp=mxp
      ios=gridrd8(gridnm,nxp,nyp,grid,x0p,x1p,y0p,y1p,z0p,z1p)
      if (ios.ne.0) stop 'addgrid: error loading grid'
      dxp=(x1p-x0p)/(nxp-1)
      dyp=(y1p-y0p)/(nyp-1)
      do iy=1,ny
	 do ix=1,nx
	    k=(iy-1)*nx+ix
	    xp=(lon0+(ix-1)*dlon-x0p)/dxp+1
	    yp=(lat0+(iy-1)*dlat-y0p)/dyp+1
	    h(k)=h(k)+gr3int(grid,xp,yp,nxp,nyp,nxp)
	 enddo
      enddo
      end
