      program geoint
      implicit none

      include "geograph.inc"
      integer l,m,p,numb(-ndeg:ndeg,0:ndeg)
      integer i,j,k,ip,rec(4)
      real*8 amp(2,-ndeg:ndeg,0:ndeg)
*     real*8 psi0(-ndeg:ndeg,0:ndeg)
      real*8 lat,lon,x0,x1,y0,y1,dx,dy,period
      real*8 x,y,z,og,og0,wm,wm0,sinwm,sininc,cosinc,t,t0,t1,dt,epoch,
     |	psi,dr
      integer cycles,nstep
      character*80 filenm
      logical xgf/.false./
      namelist /geograph_nml/ incl,period,cycles,deadband,lmax,test,
     | x0,x1,y0,y1,dx,dy
      namelist /geoint_nml/ t0,t1,dt,og0,wm0,epoch

      call getarg(1,filenm)
      call gravrd(1.0,filenm)
      call getarg(2,filenm)
      call gravrd(-1.0,filenm)

      filenm='/user/altim'
      call checkenv('ALTIM',filenm,l)
      filenm(l+1:)='/nml/geograph.nml'
      open (7,file=filenm,status='old')
      read (7,geograph_nml)
      close (7)
      filenm(l+1:)='/nml/geoint.nml'
      open (7,file=filenm,status='old')
      read (7,geoint_nml)
      close (7)
      open (7,file="geograph.nml",status="old",err=2)
      read (7,geograph_nml)
      close (7)
2     continue
      open (7,file="geoint.nml",status="old",err=3)
      read (7,geoint_nml)
      close (7)
3     continue

      nstep=nint((t1-t0)/dt)
      incl=incl*rad
      cosinc=cos(incl)
      sininc=sin(incl)
      og0=og0*rad
      wm0=wm0*rad
      n0=(2*pi)/(period*86400d0/cycles)
      wmdot=n0
      ogdot=-nint(period)*n0/cycles
      a0=(gm/n0**2)**(1d0/3d0)
      
      call getarg(3,filenm)
      if (filenm.ne.' ') then
	 open (30,file=filenm,status='new',form='unformatted',
     |   access='direct',recl=18)
	 xgf=.true.
      endif
      call d_lmp

      do ip=1,ipmax
	 l=ideg(ip)
	 m=iord(ip)
	 j=ics(ip)
	 if (m/2*2.eq.m) then
	    do p=0,l
	       k=l-2*p
	       i=pnt(l)+m*(l+1)+p
*	       psi0(k,m)=k*wmdot+m*ogdot
	       amp(j,k,m)=amp(j,k,m)+cs(ip)*dlmp(i)
	       numb(k,m)=numb(k,m)+1
	    enddo
	 else
	    do p=0,l
	       k=l-2*p
	       i=pnt(l)+m*(l+1)+p
*	       psi0(k,m)=k*wmdot+m*ogdot
	       amp(3-j,k,m)=amp(3-j,k,m)-(-1)**j*cs(ip)*dlmp(i)
	       numb(k,m)=numb(k,m)+1
	    enddo
	 endif
      enddo

      do j=0,nstep
	 t=t0+j*dt
	 wm=wm0+(t-epoch)*wmdot
	 og=og0+(t-epoch)*ogdot
	 x=cos(wm)
	 sinwm=sin(wm)
	 y=sinwm*cosinc
	 z=sinwm*sininc
	 lat=asin(z)/rad
	 lon=(atan2(y,x)+og)/rad
	 lon=lon-nint(lon/360)*360
	 dr=0
         do m=0,lmax
            do k=-lmax,lmax
	       psi=k*wm+m*og
	       dr=dr+amp(1,k,m)*cos(psi)+amp(2,k,m)*sin(psi)
	    enddo
         enddo
	 if (xgf) then
	    rec(1)=nint(t)
	    rec(2)=nint(lat*1d6)
	    rec(3)=nint(lon*1d6)
	    rec(4)=nint(dr*1d6)
	    write (30,rec=j+2) rec
	 else
	    write (*,610) t,wm/rad,og/rad,lat,lon,dr
	 endif
      enddo
      if (xgf) then
	 write (30,rec=1) '@XGF',nstep+1
      endif
  610 format (6f12.3)
      end
