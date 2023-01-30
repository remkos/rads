      subroutine gridrms(nx,ny,mx,z,n,zmin,zmax,zmean,zrms)
      real*8 z(mx,*),zmean,zrms,zmin,zmax,zz
      integer n,nx,ny,mx,kx,ky
      zmean=0
      zrms=0
      zmin=1d40
      zmax=-1d40
      n=0
      do ky=1,ny
	 do kx=1,nx
	    zz=z(kx,ky)
	    if (abs(zz).lt.1d20) then
	       zmean=zmean+zz
	       zrms=zrms+zz**2
	       zmin=min(zz,zmin)
	       zmax=max(zz,zmax)
	       n=n+1
	    endif
	 enddo
      enddo
      zmean=zmean/n
      zrms=sqrt(zrms/n)
      end
