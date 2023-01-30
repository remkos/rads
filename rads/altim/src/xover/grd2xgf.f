      program grd2xgf

      implicit none

      character*80 filenm
      integer it,ix,iy,iz,kx,ky,nx,ny,mx,irec,k
      parameter (mx=1000000)
      real grid(mx),x0,x1,y0,y1,z0,z1,dx,dy
      integer*2 is

      call getarg(1,filenm)

      nx=0
      ny=mx
      call gridrd4(filenm,nx,ny,grid,x0,x1,y0,y1,z0,z1)
      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)

      call getarg(2,filenm)
      open (10,file=filenm,status='new',access='direct',
     .form='unformatted',recl=18)

      k=0
      irec=0
      do ky=1,ny
	 iy=nint((y0+dy*(ky-1))*1e6)
	 do kx=1,nx
	    k=k+1
	    if (grid(k).lt.1e20) then
	       ix=nint((x0+dx*(kx-1))*1e6)
	       iz=nint(grid(k)*1e6)
	       is=1
	       irec=irec+1
	       write (10,rec=irec+1) it,iy,ix,iz,is
	    endif
	 enddo
      enddo
      write (10,rec=1) '@XGF',irec
      end
