      program varsgridb
************************************************************************
* Varsgridb: Gridding routine using 'distance weighting' with variable *
* sigma, incl. reinterpolation.                                        *
*                                                                      *
* The weight function, which may be changed using the statement        *
* function, equals exp(-dist**2/sigma**2), where the sigma depends on  *
* the distance of the closest measurement.			       *
*                                                                      *
* Remko Scharroo, Delft, 24 May 1991.                                  *
* 2 Dec 94 - New GRIDWR8 routine implemented.
************************************************************************
      implicit real*8 (a-h,o-z), integer (i-n)
      parameter    (MaxGridPoints=1000000)
      parameter    (MaxMeasurements=1000000)
      real*8       sum(MaxGridPoints),wgt(MaxGridPoints),
     .             dr2min(MaxGridPoints),sig2(MaxGridPoints)
      integer*4    ssh(MaxMeasurements),distyp,gridwr8
      integer*2    isg
      character*18 xgf,xgfout(MaxMeasurements)
      character    fspec*4,xgfin*80,gridut*80,text*21
      logical      reint,latcor
      equivalence  (xgf(1:4),itg),(xgf(5:8),iyg),(xgf(9:12),ixg),
     &             (xgf(13:16),izg),(xgf(17:18),isg)
*
* Statement function for weight
*
      wfunc(x2)=dexp(-x2)
*
      pi=4*atan(1d0)
      rad=pi/180
*
* Initialization for gridding steps
*
      zmin=+1d40
      zmax=-1d40
      sumz=0
      sumz2=0
*
* Read parameters and standard input
*
      call getarg(1,xgfin)
      call getarg(2,gridut)
      if (xgfin.eq.' ' .or. gridut.eq.' ') goto 1360
      read (5,*,end=1300) x0,x1,nxc
      read (5,*,end=1300) y0,y1,nyc
      read (5,*,end=1300) distyp
      read (5,*,end=1300) horfac,ulthor
      read (5,*,end=1300) wpar1,wpar2
      read (5,*,end=1300) bad
      aedit=0
      read (5,*,end=10) aedit
   10 reint=(aedit.gt.0)
      latcor=(distyp.eq.1)
      text='gridsp. (rectangular)'
      if (latcor) text='degrees (spherical)'
      write (6,1000) xgfin,gridut
*
* Values all read in. Open input XGF file.
*
      open (10,file=xgfin,status='old',form='unformatted',
     .      recl=18,access='direct',err=1307)
      read (10,rec=1) fspec,nread
      if (fspec.ne.'@XGF') goto 1350
      if (nread.gt.MaxMeasurements) goto 1308
      write (6,1010) x0,x1,nxc,y0,y1,nyc,wpar1,wpar2,text,horfac,ulthor
      if (reint) write (6,1015) aedit
*
* Perform some preparatory computations and tests.
*
      nxp=nxc+1
      nyp=nyc+1
      npt=nxp*nyp
      if (nxc.lt.1 .or. nyc.lt.1 .or. npt.gt.MaxGridPoints) goto 1302
      if (horfac.lt.1 .or. ulthor.le.0) goto 1313
      if (wpar1.le.0 .or. wpar2.lt.0 .or. wpar2.ge.10) goto 1312
      perc=100./npt
      step=(nread-.5)*.02
*
* Normalize some units by ultimate horizon
*
      fx=nxc/(x1-x0)
      fy=nyc/(y1-y0)
      bkx=ulthor
      bky=ulthor
      if (latcor) bky=bky*fy
      wpar1=wpar1**2
      wpar2=(wpar2/ulthor)**2
      horfac=horfac**2
*
* Prepare arrays for grid
*
      do 100 kpt=1,npt
  100    dr2min(kpt)=1
      ill=0
      invalid=0
      nused=0
      stip=step
      write (6,1020) nread
*
* Loop reading from file, calculate closest measurement for each grid point.
* (Pass 1). Store measurement in array xgfout().
*
      do 210 iread=1,nread
         read (10,rec=iread+1) xgf
         if (iread.ge.stip) then
            write (6,1030)
            stip=stip+step
         endif
	 if (isg.lt.0) goto 210
*
         x=ixg*1d-6
         y=iyg*1d-6
         z=izg*1d-6
*
         nused=nused+1
	 xgfout(nused)=xgf
	 if (latcor) bkx=ulthor*fx/cos(y*rad)
*
* Compute test indices for this data point
*
	 xg=(x-x0)*fx+1
	 yg=(y-y0)*fy+1
         kxlo=max(  1,int(xg-bkx+1))
         kylo=max(  1,int(yg-bky+1))
         kxhi=min(nxp,int(xg+bkx  ))
         kyhi=min(nyp,int(yg+bky  ))
*
* Check if (distance/ultimate horizon)**2 is shortest for all pertinent points
*
         do 200 ky=kylo,kyhi
            dy=(yg-ky)/bky
            do 200 kx=kxlo,kxhi
               dx=(xg-kx)/bkx
               dr2=dx**2+dy**2
               if (dr2.lt.1) then
                  kpt=kx+(ky-1)*nxp
                  dr2min(kpt)=min(dr2min(kpt),dr2)
               endif
  200    continue
  210 continue
      if (nused.eq.0) goto 1309
*
* Now compute the normalized sigma**2 for every grid point
*
      do 280 kpt=1,npt
  280    sig2(kpt)=max(wpar2,wpar1*dr2min(kpt))
*
* Loop: weights and grid values calculation (Pass 2)
*
      step=(nused-.5)*.02
      stip=step
      write (6,1025)
      do 300 iused=1,nused
         xgf=xgfout(iused)
         if (iused.ge.stip) then
            write (6,1030)
            stip=stip+step
         endif
	 if (isg.lt.0) goto 300
*
         x=ixg*1d-6
         y=iyg*1d-6
         z=izg*1d-6
*
         zmax=max(zmax,z)
         zmin=min(zmin,z)
         sumz=sumz+z
         sumz2=sumz2+z**2
*
	 if (latcor) bkx=ulthor*fx/cos(y*rad)
*
* Compute test indices for this data point
*
	 xg=(x-x0)*fx+1
	 yg=(y-y0)*fy+1
         kxlo=max(  1,int(xg-bkx+1))
         kylo=max(  1,int(yg-bky+1))
         kxhi=min(nxp,int(xg+bkx  ))
         kyhi=min(nyp,int(yg+bky  ))
*
         do 290 ky=kylo,kyhi
            dy=(yg-ky)/bky
            do 290 kx=kxlo,kxhi
               dx=(xg-kx)/bkx
               kpt=kx+(ky-1)*nxp
               hori=horfac*sig2(kpt)
               dr2=dx**2+dy**2
               if (dr2.le.hori) then
                  w=wfunc(dr2/sig2(kpt))
                  wgt(kpt)=wgt(kpt)+w
                  sum(kpt)=sum(kpt)+w*z
               endif
  290    continue
  300 continue
*
* Postprocess grid
*
      zmean=sumz/nused
      zsigma=dsqrt((sumz2-sumz**2/nused)/(nused-1))
*
      do 350 kpt=1,npt
         w=wgt(kpt)
         if (w.lt..1) ill=ill+1
         if (w.ne.0) then
            sum(kpt)=sum(kpt)/w
         else
            sum(kpt)=bad
            invalid=invalid+1
         endif
  350 continue
*
* Now start output of grid
*
      if (gridwr8(gridut,nxp,nyp,sum,nxp,x0,x1,y0,y1).ne.0)
     &          write (6,550) 'grid probably unusable'
      if (.not.reint) then
         write (6,1070) nread,zmin,zmean,nused,zmax,zsigma,
     .   npt,ill,ill*perc,invalid,invalid*perc
         goto 9999
      endif
*
* Initialization for re-interpolation
*
      rmin=+1d40
      rmax=-1d40
      sumz=0
      sumz2=0
      nwrit=0
*
* Reinterpolation pass 1
*
      stip=step
      write (6,1026)
      do 400 iused=1,nused
         xgf=xgfout(iused)
         if (iused.ge.stip) then
            write (6,1030)
            stip=stip+step
         endif
*
         x=ixg*1d-6
         y=iyg*1d-6
         z=izg*1d-6
*
         xg=(x-x0)*fx
         yg=(y-y0)*fy
         resid=z-gr1int(sum,xg,yg,nxp,nyp,nxp)
         ssh(iused)=nint(resid*1d6)
*
* Keep statistics
*
	 rmax=max(rmax,resid)
	 rmin=min(rmin,resid)
         sumz=sumz+resid
         sumz2=sumz2+resid**2
  400 continue
*
* Prepare edit criteria
*
      rmean=sumz/nused
      rsigma=dsqrt((sumz2-sumz**2/nused)/(nused-1))
      rlo=rmean-aedit*rsigma
      rhi=rmean+aedit*rsigma
      minedt=nint(rlo*1d6)
      maxedt=nint(rhi*1d6)
*
* Edit and output. Re-open and rewrite XGF-file. Reinterpolation pass 2.
*
      close (10,status='delete')
      open (10,file=xgfin,status='old',form='unformatted',
     .      recl=18,access='direct',err=1307)
      stip=step
      write (6,1027)
      do 500 iused=1,nused
	 xgf=xgfout(iused)
         if (iused.ge.stip) then
            write (6,1030)
            stip=stip+step
         endif
         izg=ssh(iused)
         if (izg.ge.minedt .and. izg.le.maxedt) then
            nwrit=nwrit+1
            write (10,rec=nwrit+1) xgf
         endif
  500 continue
*
      write (10,rec=1) '@XGF',nwrit
      close (10)
      nedit=nused-nwrit
      write (6,1070) nread,zmin,zmean,nused,zmax,zsigma,
     .npt,ill,ill*perc,invalid,invalid*perc
      write (6,1080) rmin,rmean,rmax,rsigma,nedit,rlo,rhi,nwrit
      goto 9999
*
* Formats
*
  550 format (a)
 1000 format (79('*')/'*',21x,'variable distance-weighted gridding',
     .t79,'*'/79('*')/' XGF : ',a70/'Grid : ',a70/)
 1010 format ('x-interval: ',2f7.1,4x,'nr. of cells:',i4/
     .'y-interval: ',2f7.1,4x,'nr. of cells:',i4//
     .'Gridding:   sigma=',f7.3,'*dist  with a minimum of',f7.3,1x,a21/
     .10x,'horizon=',f7.3,'*sigma with a maximum of',f7.3/)
 1015 format ('Reinterpolation: edit multiplier on residual sigma:',
     .f7.3/)
 1020 format(17x,
     .'2%  .    .    .    .   50%   .    .    .    .  100% (',
     .i8,')'/' Gridding Pass 1 ',$)
 1025 format(/' Gridding Pass 2 ',$)
 1026 format (/'Reinterp. Pass 1 ',$)
 1027 format (/'Reinterp. Pass 2 ',$)
 1030 format ('#',$)
 1070 format (//i8,' data points read       min:',f9.3,5x,'mean:',f9.3/
     .i8,' data points used       max:',f9.3,4x,'sigma:',f9.3//
     .i8,' grid points',i10,' ill-determined (',f4.1,'%)'/
     .i30,'  undetermined  (',f4.1,'%)'/)
 1080 format('Residuals: min:',f9.3,5x,'mean:',f9.3/
     .11x,'max:',f9.3,4x,'sigma:',f9.3//
     .i8,' data points edited     min:',f9.3,6x,'max:',f9.3/
     .i8,' data points written')
*
* Errors
*
 1300 write (6,550) 'premature EOD on standard input.'
      goto 9999
*
 1302 write (6,550) 'illegal grid dimensions.'
      goto 9999
*
 1307 write (6,550) 'error opening XGF-file.'
      goto 9999
*
 1308 write (6,550) 'too many data points in XGF file.'
      goto 9999
*
 1309 write (6,550) 'no data to process.'
      goto 9999
*
 1312 write (6,550) 'illegal weight factor or low weight boundary.'
      goto 9999
*
 1313 write (6,550) 'horizon must be > 1*sigma and > 0 degrees.'
      goto 9999
*
 1350 write (6,550) 'input data not in XGF format.'
      goto 9999
*
 1360 write (6,550) 'usage: varsgridb <XGF file> <grid file>'
 9999 end
