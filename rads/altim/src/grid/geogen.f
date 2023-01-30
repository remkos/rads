c
c     compute geoid heights or gravity anomalies
c
c      ichoice:
c      0--> geocentric
c      1--> geodetic  
c      2--> geodetic on geocentric grid
c
c   
      program geogen

      implicit real*8 (a-h,o-z), integer (i-n)
      integer maxdeg,nlab,nphi,nmax,number
      parameter(maxdeg=360)
      parameter(nlab=3621,nphi=1801)
      parameter(nmax=360)
      parameter(number=(nmax+1)*(nmax+2)/2)
      
      character*80 input,output,arg
      
      logical new,ano

      real*8 pnm(0:maxdeg,0:maxdeg)
      real*8 cnm(number),snm(number),cap(0:6)
      real*8 am(0:maxdeg),bm(0:maxdeg),rath

      real grid(nlab,nphi),rlonmin,rlonmax,rlatmin,rlatmax

      integer offset(0:nmax)
      
      call getarg(0,arg)
      new=index(arg,'geogen-new').gt.0
      call getarg(1,input)
      call getarg(2,output)
      if (output.eq.' ') then
	 write (0,600)
	 goto 9999
      endif
600   format('usage: geogen-new <gravity model> <grid file>'/
     .'   or: geogen <gravity model> <grid file> [-hgt | -ano]')

      if (new) then
	 write(6,551) 'geoid heights (0) or gravity anomalies (1) --> '
	 read(5,*) ichoise
	 ano=(ichoise.eq.1)
      else
	 call getarg(3,arg)
	 ano=(arg.eq.'-ano')
      endif

      write(6,'(a,a,$)')
     .'geocentric(0), geodetic(1),',
     .' or geodetic values on geocentric grid (2) -> '
      read(5,*) ichoice
       
      write(6,551) 'minimum and maximum longitude (degrees) -> '
      read(5,*) rlonmin,rlonmax
      
      write(6,551) 'minimum and maximum latitude (degrees) -> '
      read(5,*) rlatmin,rlatmax
      
      write(6,551) 'number of cells per degree (lon/lat) -> '
      read(5,*) rnrlon,rnrlat
      steplon=1d0/rnrlon
      steplat=1d0/rnrlat
   
      nlab2=1+nint((rlonmax-rlonmin)/steplon)
      nphi2=1+nint((rlatmax-rlatmin)/steplat)

      if (nlab2.gt.nlab .or. nphi2.gt.nphi) then
	 write (6,*) 'ERROR: grid too large',nlab2,nphi2
	 stop
      endif
      
      write(6,551) 'minimum and maximum degree of gravity field -> '
      read(5,*) lmin,lmax
c      if (lmin.lt.2) lmin=2
      lmin=max(1,lmin)
      
  551 format(a,$)

      open(unit=10,status='old',form='unformatted',file=input)

* Read osu91a model

      read(10) (cnm(iadr),iadr=1,number)
      read(10) (snm(iadr),iadr=1,number)
    
* Reference ellipsoid

      ae=6378137.0d0
      r=ae
      f=1d0/298.257d0
      f2=(1d0-f)**2
      fmin=2d0*f-f**2
      gm=3.986004404d14
      cap(0)=1d0
      cap(2)=-4.8416542d-04
      cap(4)=7.90297553d-07
      cap(6)=-1.6872191d-09

      do m=0,nmax
	 offset(m) = m*(nmax+1) - m*(m+1)/2 + 1
      enddo
  
* Correct coefficients for reference ellipsoid

      cnm(offset(0)+0)=cnm(offset(0)+0)-cap(0)
      cnm(offset(0)+2)=cnm(offset(0)+2)-cap(2)
      cnm(offset(0)+4)=cnm(offset(0)+4)-cap(4)
      cnm(offset(0)+6)=cnm(offset(0)+6)-cap(6)
      
      pi=4*datan(1d0)
      rad=pi/180
      rnobs=0d0
      rmn=0d0
      rms=0d0
      nobs=nlab*nphi

* Compute geoid heights or gravity anomalies

      call statbar(0,nlab2*nphi2,'Generating grid')
      do iphi=1,nphi2
	 phi=(rlatmin+(iphi-1)*steplat)*rad
	 if (ichoice.eq.0) then
	    r=ae
	 else if (ichoice.eq.1) then
	    phi=datan(dtan(phi)*f2)
	    rhelp=ae**2*f2/(1d0-fmin*dcos(phi)**2)
	    r=dsqrt(rhelp)
	 else if (ichoice.eq.2) then
	    rhelp=ae**2*f2/(1d0-fmin*dcos(phi)**2)
	    r=dsqrt(rhelp)
	 endif
	 gamma=gm/r**2
	 ratio=ae/r
	 ct=dcos(pi/2-phi)
	 st=dsin(pi/2-phi)
	 call legpol(ct,st,pnm,lmax)
	 do m=0,lmax
	    am(m)=0.d0
	    bm(m)=0.d0
	 enddo
	 if (ano) then
	    factor=gamma*1d5
	 else
	    factor=r
	 endif
	 do l=lmin,lmax
	    rath=ratio**l
	    do m=0,l
	       iadr=offset(m)+l
	       if (ano) then
	          am(m)=am(m)+(l-1)*pnm(l,m)*cnm(iadr)*rath
	          bm(m)=bm(m)+(l-1)*pnm(l,m)*snm(iadr)*rath
	       else
	          am(m)=am(m)+pnm(l,m)*cnm(iadr)*rath
	          bm(m)=bm(m)+pnm(l,m)*snm(iadr)*rath
	       endif
	    enddo
	 enddo
	 do ilab=1,nlab2
	    rlon=(rlonmin+(ilab-1)*steplon)*rad
	    value=0.0d0
	    do m=0,lmax
	       value=value+am(m)*dcos(m*rlon)+bm(m)*dsin(m*rlon)
	    enddo
	    value=value*factor
	    grid(ilab,iphi)=value
	    rms=rms+(st*value)**2
	    rmn=rmn+st*value
	    rnobs=rnobs+st
	 enddo
	 call statbar(1,iphi*nlab2,' ')
      enddo

* Write geoid heights or gravity anomalies

      call gridwr4(output,nlab2,nphi2,grid,nlab,
     &  rlonmin,rlonmax,rlatmin,rlatmax)
      
* Write statistics

      if (ano) then
	 write(6,610) rmn/rnobs,dsqrt(rms/rnobs)
      else
	 write(6,620) rmn/rnobs,dsqrt(rms/rnobs)
      endif


  610 format ('Mean anomaly (mgal): ',f10.3/
     .        ' RMS anomaly (mgal): ',f10.3)
  620 format ('Mean geoid height (m): ',f10.3/
     .        ' RMS geoid height (m): ',f10.3)

9999  end
