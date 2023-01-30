c
c     compute geoid heights or gravity anomalies
c   
      program geogen3

      implicit none
      include "geograph.inc"
      integer*4 maxdeg,mx,my,nmax,number
      parameter(maxdeg=360,nmax=maxdeg)
      parameter(mx=3621,my=1801)
      parameter(number=(nmax+1)*(nmax+2)/2)
      
      integer*4 lmin/1/,iadr,i,kx,ky,nx,ny,ios,l,ichoice/1/,m,iargc
      character*80 input1/'-'/,input2/'-'/,output/'-'/,arg
      
      logical ano/.false./

      real*8 cnm(number),snm(number),cap(0:6)
      real*8 am(0:maxdeg),bm(0:maxdeg),rath,grid(mx,my),plm(0:maxdeg)
      real*8 lon0/-180d0/,lon1/180d0/,dlon/1d0/
      real*8 lat0/-90d0/,lat1/90d0/,dlat/1d0/
      real*8 f,f2,fmin,r/0d0/,coslat,ratio,value,wmean,wrms,wnr,
     |       lon,lat,gamma,factor

      integer offset(0:nmax)

      ae=6378137.0d0
      f=1d0/298.257d0
      f2=(1d0-f)**2
      fmin=2d0*f-f**2
      gm=3.986004404d14
      lmax=70
      
      do i=1,iargc()
         call getarg(i,arg)
	 if (arg(:1).eq.'-a') then
	    ano=.true.
	 else if (arg(:4).eq.'geo=') then
	    read (arg(5:),*) ichoice
	 else if (arg(:4).eq.'alt=') then
	    read (arg(5:),*) r
	    ichoice=3
	 else if (arg(:4).eq.'deg=') then
	    read (arg(5:),*,iostat=ios) lmin,lmax
	    lmin=max(1,lmin)
	 else if (arg(:4).eq.'lon=') then
	    read (arg(5:),*,iostat=ios) lon0,lon1,dlon
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*,iostat=ios) lat0,lat1,dlat
	 else if (arg(:5).eq.'grid=') then
	    output=arg(6:)
	 else if (input1.eq.'-') then
	    input1=arg
	 else
	    input2=arg
	 endif
      enddo

      if (input1.eq.'-') then
	 write (0,600)
	 goto 9999
      endif
600   format('geogen3 -- Generate geoid (difference) grid'//
     |'usage: geogen3 [options] model1 [model2]'//
     |'where:'/
     |' model1       : gravity field model file'/
     |' model2       : optional gravity field model file',
     |' (for difference)'//
     |'and [options] are:'/
     |' -ano         : produce anomalies in stead of geoid heights'/
     |' geo=i        : i specifies reference frame (0=geocentric,',
     |' 1=geodetic, 2=geodetic on'/
     |'                geocentric grid, 3=at altitude) (default: 1)'/
     |' alt=h        : specifies altitude in km (implies geo=3)'/
     |' deg=low,hi   : specify minimum and maximum degree of gravity',
     |' field (default: 1,70)'/
     |' grid=name    : specify name of output grid'/
     |' lon=x0,x1,dx : longitude boundaries and spacing',
     |' (default: -180,180,1)'/
     |' lat=y0,y1,dy : latitude boundaries and spacing',
     |' (default: -90,90,1)')

      nx=1+nint((lon1-lon0)/dlon)
      ny=1+nint((lat1-lat0)/dlat)

      if (nx.gt.mx .or. ny.gt.my) then
	 write (*,*) 'ERROR: grid too large',nx,ny
	 stop
      endif
      
      call gravrd(1.0,input1)
      if (input2.ne.'-') call gravrd(-1.0,input2)

      do m=0,nmax
	 offset(m) = m*(nmax+1) - m*(m+1)/2 + 1
      enddo
      do i=1,ipmax
         iadr=offset(iord(i))+ideg(i)
	 if (ics(i).eq.1) then
	    cnm(iadr)=cs(i)
	 else
	    snm(iadr)=cs(i)
	 endif
      enddo
    
* Reference ellipsoid

      if (input2.eq.'-') then
         cap(0)=1d0
         cap(2)=-4.8416542d-04
         cap(4)=7.90297553d-07
         cap(6)=-1.6872191d-09

         cnm(offset(0)+0)=cnm(offset(0)+0)-cap(0)
         cnm(offset(0)+2)=cnm(offset(0)+2)-cap(2)
         cnm(offset(0)+4)=cnm(offset(0)+4)-cap(4)
         cnm(offset(0)+6)=cnm(offset(0)+6)-cap(6)
      endif
      
      wnr=0d0
      wmean=0d0
      wrms=0d0

* Compute geoid heights or gravity anomalies

      call statbar(0,nx*ny,'Generating grid')
      do ky=1,ny
	 lat=(lat0+(ky-1)*dlat)*rad
	 if (ichoice.eq.0) then
	    r=ae
	 else if (ichoice.eq.1) then
	    lat=atan(tan(lat)*f2)
	    r=ae**2*f2/(1d0-fmin*cos(lat)**2)
	    r=sqrt(r)
	 else if (ichoice.eq.2) then
	    r=ae**2*f2/(1d0-fmin*cos(lat)**2)
	    r=sqrt(r)
	 else
	    r=r*1d3+ae
	 endif
	 gamma=gm/r**2
	 ratio=ae/r
	 coslat=cos(lat)
	 do l=0,lmax
	    am(l)=0.d0
	    bm(l)=0.d0
	 enddo
	 if (ano) then
	    factor=gamma*1d5
	 else
	    factor=r
	 endif
	 do l=lmin,lmax
	    rath=ratio**l
	    call p_lm(l,lat,plm)
	    do m=0,l
	       iadr=offset(m)+l
	       if (ano) then
	          am(m)=am(m)+(l-1)*plm(m)*cnm(iadr)*rath
	          bm(m)=bm(m)+(l-1)*plm(m)*snm(iadr)*rath
	       else
	          am(m)=am(m)+plm(m)*cnm(iadr)*rath
	          bm(m)=bm(m)+plm(m)*snm(iadr)*rath
	       endif
	    enddo
	 enddo
	 do kx=1,nx
	    lon=(lon0+(kx-1)*dlon)*rad
	    value=0.0d0
	    do m=0,lmax
	       value=value+am(m)*cos(m*lon)+bm(m)*sin(m*lon)
	    enddo
	    value=value*factor
	    grid(kx,ky)=value
	    wrms=wrms+coslat*value**2
	    wmean=wmean+coslat*value
	    wnr=wnr+coslat
	 enddo
	 call statbar(1,ky*nx,' ')
      enddo

* Write geoid heights or gravity anomalies

      if (output.ne.'-') then
         call gridwr8(output,nx,ny,grid,mx,
     |     lon0,lon1,lat0,lat1)
      endif
      
* Write statistics

      if (ano) then
	 write(*,610) wmean/wnr,sqrt(wrms/wnr)
      else
	 write(*,620) wmean/wnr,sqrt(wrms/wnr)
      endif


  610 format ('Mean anomaly (mgal): ',f10.3/
     |        ' RMS anomaly (mgal): ',f10.3)
  620 format ('Mean geoid height (m): ',f10.3/
     |        ' RMS geoid height (m): ',f10.3)

9999  end
