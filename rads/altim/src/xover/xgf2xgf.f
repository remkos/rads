      program xgf2xgf
*
* Program xgf2xgf converts XGF files to XGF or makes changes
* to XGF files.
*
*-
*  3-Nov-1999 - Added DATEARG function
* 28-Jul-2000 - Added fit= flag
*  4-Jun-2003 - Better initialisation of variables
*----------------------------------------------------------------------
      implicit none
      integer*4 type,mode/0/,field/0/,mrec,nrec,irec,
     |          narg,iargc,iarg,if,nf,nx,ny,mx,istep
      integer*4 trk0,trk1,it0,it1,l0,l1
      parameter (mx=3601*1801)
      integer*4 xgf4(4),xaf4(7),xxf4(12)
      integer*2 xgf2(9),xaf2(14),xxf2(24)
      real*8    grid(mx),lat,lon,x,y,dx/1d0/,dy/1d0/,
     |		gfact/-1d6/,hfact/1d0/, gr3int,glon0,glon1,glat0,glat1,
     |		rmin,rmax,rmean,rrms,r,h,nfact,
     |		lon0/-180d0/,lon1/180d0/,lat0/-90d0/,lat1/90d0/,
     |		t,t0,t1,dum,fit(4)/4*0d0/,rad
      real*4	noise
      integer*4 ios,openf,readf,writef,closef,seekf,fdin,fdout
      logical   ref/.false./,addnoise/.false./,sum/.false./,datearg,
     |		usefit/.false./
      integer*4 minint4,maxint4,mf
      parameter (maxint4=2147483647,minint4=-maxint4-1,mf=400)
      character*80 filenm(mf),arg,spec*4,gridnm
      equivalence (xgf4,xgf2),(xgf4,xaf4),(xgf4,xaf2),(xxf4,xxf2)

* Initialize

      istep=1
      nf=0
      it0=minint4
      it1=maxint4
      trk0=minint4
      trk1=maxint4
      narg=iargc()
      rad=atan(1d0)/45

      call getarg(0,arg)
      if (index(arg,'xaf2xgf').ne.0) mode=1
      if (index(arg,'xxf2xgf').ne.0) mode=2

      if (narg.eq.0) then
	 if (mode.eq.0) write (0,620)
	 if (mode.eq.1) write (0,621)
	 if (mode.eq.2) write (0,622)
	 write (0,623)
	 if (mode.ge.1) write (0,624)
	 if (mode.eq.2) write (0,625)
	 goto 9999
      endif

      do iarg=1,narg
         call getarg(iarg,arg)
	 if (arg(:5).eq.'step=') then
	    read (arg(6:),*) istep
	    write (0,*) 'Step = ',istep
	 else if (arg(:4).eq.'ref=') then
	    gridnm=arg(5:)
	    nx=0
	    ny=mx
	    call gridrd8(gridnm,nx,ny,grid,glon0,glon1,glat0,glat1,
     |		rmin,rmax)
	    dx=(glon1-glon0)/(nx-1)
	    dy=(glat1-glat0)/(ny-1)
	    ref=.true.
         else if (arg(:6).eq.'gfact=') then
            read (arg(7:),*) gfact
	    gfact=gfact*1d6
         else if (arg(:6).eq.'hfact=') then
            read (arg(7:),*) hfact
         else if (arg(:6).eq.'noise=') then
            read (arg(7:),*) nfact
	    nfact=nfact*1d6
	    addnoise=.true.
         else if (arg(:4).eq.'fit=') then
            read (arg(5:),*) fit
            usefit=.true.
	 else if (arg(:4).eq.'lon=') then
	    read (arg(5:),*) lon0,lon1
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*) lat0,lat1
	 else if (datearg(arg,t0,t1,dum)) then
	    it0=nint(t0)
	    it1=nint(t1)
	 else if (arg(:4).eq.'trk=') then
	    read (arg(5:),*) trk0,trk1
	 else if (arg(:2).eq.'-a') then
	    field=0
	 else if (arg(:2).eq.'-p') then
	    field=1
	 else if (arg(:2).eq.'-d') then
	    field=2
	 else if (arg(:2).eq.'-m') then
	    sum=.true.
	 else if (arg(:2).eq.'-x') then
	    sum=.false.
	 else
	    nf=nf+1
	    if (nf.gt.mf) call fin('too many input files')
	    filenm(nf)=arg
	 endif
      enddo

      fdout=openf(filenm(nf),'w')
      if (fdout.le.0) call fin('xgf2xgf: error opening output file')
      ios=writef(fdout,18,xgf2)
      nf=nf-1
      mrec=0
      rmean=0
      rrms=0
      rmin=1d40
      rmax=-1d40
      do if=1,nf
	 fdin=openf(filenm(if),'r')
	 ios=readf(fdin,4,spec)
	 ios=readf(fdin,4,nrec)
         l0=1
         l1=nrec
	 if (spec.eq.'@XGF') then
	    type=1
	    ios=seekf(fdin,18,0)
	 else if (spec.eq.'@XAB') then
	    type=2
	    call xtflimits(filenm(if),it0,it1,trk0,trk1,l0,l1)
	    ios=seekf(fdin,28*l0,0)
	 else if (spec.eq.'@XXB') then
	    type=3
	    ios=seekf(fdin,48,0)
	 else if (spec.eq.'@XXO') then
	    type=4
	    ios=seekf(fdin,44,0)
	 else
	    write (0,680)
	    goto 9999
	 endif
	 do irec=l0,l1,istep
	    if (type.eq.1) then
	       ios=readf(fdin,18,xgf2)
	       if (xgf2(9).lt.0) goto 100
	    else if (type.eq.2) then
	       ios=readf(fdin,28,xaf2)
	       if (xaf2(13).lt.0) goto 100
	       if (xaf2(14).lt.trk0 .or. xaf2(14).gt.trk1) goto 100
	       if (field.eq.1) then
	          xgf4(4)=xaf4(5)
	       else if (field.eq.2) then
	          xgf4(4)=xaf4(4)-xaf4(5)
	       endif
	       xgf2(9)=xaf2(13)
	    else if (type.eq.3) then
	       ios=readf(fdin,48,xxf2)
	       if (xxf2(9).lt.trk0 .or. xxf2(9).gt.trk1) goto 100
	       if (xxf2(10).lt.trk0 .or. xxf2(10).gt.trk1) goto 100
	       if (field.eq.0) then
		  if (sum) then
		     xgf4(4)=(xxf4(6)+xxf4(7))/2
		  else
		     xgf4(4)=xxf4(6)-xxf4(7)
		  endif
	       else if (field.eq.1) then
	          if (xxf2(23).lt.0 .or. xxf2(24).lt.0) goto 100
		  if (sum) then
		     xgf4(4)=(xxf4(8)+xxf4(9))/2
		  else
		     xgf4(4)=xxf4(8)-xxf4(9)
		  endif
	       else
	          if (xxf2(23).lt.0 .or. xxf2(24).lt.0) goto 100
		  if (sum) then
		     xgf4(4)=(xxf4(6)-xxf4(8)+xxf4(7)-xxf4(9))/2
		  else
		     xgf4(4)=(xxf4(6)-xxf4(8))-(xxf4(7)-xxf4(9))
		  endif
               endif
	       xgf4(2)=xxf4(1)
	       xgf4(3)=xxf4(2)
	       xgf4(1)=(xxf4(3)+xxf4(4))/2
	       xgf2(9)=(xxf2(23)+xxf2(24))/2
	    else if (type.eq.4) then
	       ios=readf(fdin,44,xxf2)
	       if (sum) then
		  xgf4(4)=(xxf4(8)+xxf4(9))/2
	       else
		  xgf4(4)=xxf4(8)-xxf4(9)
	       endif
	       xgf4(2)=xxf4(1)
	       xgf4(3)=xxf4(2)
	       xgf4(1)=(xxf4(3)+xxf4(4))/2
	       xgf2(9)=0
	    endif
	    t=xgf4(1)
	    lat=xgf4(2)/1d6
	    lon=xgf4(3)/1d6
	    if (t.lt.it0 .or. t.gt.it1) goto 100
	    if (lon.lt.lon0) lon=lon+360
	    if (lon.gt.lon1) lon=lon-360
	    if (lon.lt.lon0) goto 100
	    if (lat.lt.lat0 .or. lat.gt.lat1) goto 100

	    h=hfact*xgf4(4)
	    if (ref) then
	       y=(lat-glat0)/dy+1
	       x=(lon-glon0)/dx+1
	       r=gr3int(grid,x,y,nx,ny,nx)
	       if (abs(r).gt.1d20) goto 100
	       h=h+gfact*r
	    endif
	    if (addnoise) then
	       h=h+nfact*noise()
	    endif
            if (usefit) then
               h=h-(fit(1)+
     |           cos(lat*rad)*(fit(2)*cos(lon*rad)+fit(3)*sin(lon*rad))+
     |           fit(4)*sin(lat*rad))*1d6
	    endif
	    xgf4(4)=nint(h)

	    mrec=mrec+1
	    ios=writef(fdout,18,xgf2)
	    rmin=min(rmin,h)
	    rmax=max(rmax,h)
	    rmean=rmean+h
	    rrms=rrms+h**2
100         continue
	 enddo
	 ios=closef(fdin)
      enddo

      spec='@XGF'
      ios=seekf(fdout,0,0)
      ios=writef(fdout,4,spec)
      ios=writef(fdout,4,mrec)
      ios=closef(fdout)
      rmean=rmean/mrec/1d6
      rrms=sqrt(rrms/mrec)/1d6
      write (0,600) mrec,rmin/1d6,rmax/1d6,
     |		rmean,rrms,sqrt(rrms**2-rmean**2)
      goto 9999

600   format ('Records selected :',i10/
     |'Minimum : ',f10.4/
     |'Maximum : ',f10.4/
     |'Mean    : ',f10.4/
     |'RMS     : ',f10.4/
     |'Sigma   : ',f10.4)
620   format ('xgf2xgf converts XGF files to XGF files'//
     |'usage: xgf2xgf [options] XGF(s) XGF'/)
621   format ('xaf2xgf converts XAF files to XGF files'//
     |'usage: xaf2xgf [options] XAF(s) XGF'/)
622   format ('xxf2xgf converts XXF files to XGF files'//
     |'usage: xxf2xgf [options] XXF(s) XGF'/)
623   format ('where [options] are:'/
     |'ref=gridnm      : use grid as reference'/
     |'gfact=FACTOR    : multiply reference by FACTOR (def:-1)'/
     |'hfact=FACTOR    : multiply XGF data by FACTOR (def:1)'/
     |'noise=VAR       : add Gaussian noise with variance VAR (def:0)'/
     |'fit=dR,dX,dY,dZ : remove effect of coordinate shift (m) from',
     |' sea surface heights (output of ''earthfit'')'/
     |'step=STEP       : take each out of STEP points'/
     |'lon=lon0,lon1   : select longitude boundaries'/
     |'lat=lat0,lat1   : select latitude  boundaries'/
     |'t=t0,t1         : select time frame'/
     |'              ... or use mjd=, doy=, ymd=, sec=') 
624   format(
     |'trk=trk0,trk1   : select on track number'/
     |'-a              : select a priori ssh field (default)'/
     |'-p              : select a posteriori ssh field (after inimini)'/
     |'-d              : select orbit adjustment from inimini')
625   format(
     |'-m              : select mean of asc and des'/
     |'-x              : select diff of asc and des (default)')
680   format ('xgf2xgf: input file not XGF/XXB/XAB/XXO')

9999  close (20)
      end
