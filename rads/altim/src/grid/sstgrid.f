      program sstgrid
*------------------------------------------------------------------------
*
*     SSTGRID computes a rectangular grid from SST coefficients 
*     (Legendre Polynomials). The input is a COEF file, and the output
*     is a user-defined grid
*    
*     copyright 1992/1993 by Marc Naeije and Pieter Visser
*     copyright 1995 Naeije; adapted to new grid formats, and
*                    geodetic to geocentric conversion added
*
*------------------------------------------------------------------------
      implicit none

      integer*4    ndeg,mcoef
      parameter    (ndeg=36,mcoef=(ndeg+1)*(ndeg+2)/2)
      integer*4    l,m,iadr,offset(0:ndeg)
      real*8       pnm(0:ndeg,0:ndeg),c,s,cnm(mcoef),snm(mcoef),
     |		   am(0:ndeg),bm(0:ndeg)

      character*80 arg,gridnm
      integer*4    mx,nx,ny,kx,ky,gridwr4
      parameter    (mx=1449*721)
      real*4	   grid(mx),x0,x1,y0,y1
      integer*4    k,iarg,iargc,l0,l1,lmax,ios

      real*8       rms,rmin,rmax,pi,rad,ct,st,x,y,dx,dy,height,
     |		   f,f2
      parameter    (f=1d0/298.257d0,f2=(1-f)**2)

* Initialisations

      gridnm=' '
      rms=0d0
      pi=4*datan(1d0)
      rad=pi/180d0
      rmax=-1d35
      rmin=1d35
      x0=-180
      x1=180
      dx=1
      y0=-90
      y1=90
      dy=1
      l0=0
      l1=ndeg
      lmax=0

      do m=0,ndeg
	 offset(m)=m*(ndeg+1)-m*(m+1)/2+1
      enddo

* Get user-defined input info

      do iarg=1,iargc()
	 call getarg(iarg,arg)
	 if (arg(1:4).eq.'lat=') then
	    read (arg(5:),*,iostat=ios) y0,y1,dy
	 else if (arg(1:4).eq.'lon=') then
	    read (arg(5:),*,iostat=ios) x0,x1,dx
	 else if (arg(1:2).eq.'l=') then
	    read (arg(3:),*) l0,l1
	 else
	    gridnm=arg
	 endif
      enddo

* Check if at least grid name is given

      if (gridnm.eq.' ') then
	 write (6,600)
	 stop
      endif
600   format ('sstgrid: compute SST grid from coefficients on stdin'//
     |'syntax: sstgrid [ options ] gridname'//
     |'with [ options ]'/
     |'l=l0,l1      : use only degree l0 to l1 (def: all)'/
     |'lon=x0,x1,dx : longitude boundaries and stepsize ',
     |'(def: -180,180,1)'/
     |'lat=y0,y1,dy : latitude  boundaries and stepsize ',
     |'(def: -90,90,1)')
	  
* Proces input info

      nx=nint((x1-x0)/dx+1)
      ny=nint((y1-y0)/dy+1)
      if (nx*ny.gt.mx) stop "sstgrid: grid too large"
      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)
      
* Read Input coefficient file

10    read (5,*,end=20) l,m,c,s
      if (l.gt.ndeg) stop "sstgrid: degree/order too high"
      lmax=max(lmax,l)
      iadr=offset(m)+l
      cnm(iadr)=c
      snm(iadr)=s
      goto 10
20    continue

* Fill grid

      k=0
      l1=min(l1,lmax)
      call statbar(0,nx*ny,'Generating grid')
      do ky=1,ny
         y=(y0+(ky-1)*dy)*rad
	 y=datan(dtan(y)*f2)
         ct=dcos(pi/2-y)
         st=dsin(pi/2-y)
         call legpol(ct,st,pnm,l1)
	 do m=0,l1
	    am(m)=0
	    bm(m)=0
	 enddo
	 do l=l0,l1
	    do m=0,l
	       iadr=offset(m)+l
	       am(m)=am(m)+pnm(l,m)*cnm(iadr)
	       bm(m)=bm(m)+pnm(l,m)*snm(iadr)
	    enddo
	 enddo
         do kx=1,nx
	    k=k+1
            height=0
            x=(x0+(kx-1)*dx)*rad
            do m=0,l1
               height=height+am(m)*cos(m*x)+bm(m)*sin(m*x)
	    enddo
            rms=rms+height**2
            rmax=max(rmax,height)
            rmin=min(rmin,height)
            grid(k)=height
	 enddo
	 call statbar(1,ky*nx,' ')
      enddo

      if (gridwr4(gridnm,nx,ny,grid,nx,x0,x1,y0,y1).gt.0)
     |		stop "sstgrid: error occurred while writing grid"
      write(6,302) rmin,rmax,dsqrt(rms/k)
302   format('Min (m): ',f8.3,' Max (m): ',f8.3,' RMS (m): ',f8.3)

      end
