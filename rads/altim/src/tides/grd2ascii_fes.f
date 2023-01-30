      program grd2ascii

      implicit none
      integer*4 nxmax,nymax,nx,ny,kx,ky,iargc,i
      parameter (nxmax=721,nymax=361)
      real*8 xmin,ymin,xmax,ymax,zmin,zmax,dx,dy,mask
      real*8 wra(nxmax,nymax),wrg(nxmax,nymax)
      character*80 ampname,phaname

      if (iargc().ne.2) then
	 write (*,500)
	 goto 9999
      endif
500   format('grd2ascii -- Convert DEOS files for old FES to ASCII'//
     |'syntax: grd2ascii <ampname> <phaname>')

      call getarg(1,ampname)
      call getarg(2,phaname)

* Read input parameters

      nx=nxmax
      ny=nymax
      call gridrd8(ampname,nx,ny,wra,xmin,xmax,ymin,ymax,zmin,zmax)
      call gridrd8(phaname,nx,ny,wrg,xmin,xmax,ymin,ymax,zmin,zmax)
      dx=(xmax-xmin)/(nx-1)
      dy=(ymax-ymin)/(ny-1)
      mask=999.9

      write (*,600) xmin,xmax,ymin,ymax,dx,dy,nx,ny,mask,mask
600   format(2f9.3/2f9.3/2f9.3/2i6/2f9.3)

      do ky=1,ny
         do kx=1,nx
	    if (wra(kx,ky).ge.mask .or. wrg(kx,ky).ge.mask) then
	       wra(kx,ky)=mask
	       wrg(kx,ky)=mask
	    endif
	 enddo
	 do kx=1,nx,30
	    write (*,'(30f7.2)') (wra(i,ky),i=kx,min(kx+29,nx))
	    write (*,'(30f7.1)') (wrg(i,ky),i=kx,min(kx+29,nx))
	 enddo
      enddo
9999  end
