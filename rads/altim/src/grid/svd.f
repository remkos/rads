      implicit none
      integer nx,ny,mmax,nmax,nobs,nargs,iargc,i,j,k,
     .        neof,lutgrid,lnblnk,maskx,masky,amax,ndim
      parameter (mmax=361*181,nmax=200,amax=mmax*150)
      real*8 a(amax),v(nmax*nmax),w(nmax),grid(mmax),
     .       mask(mmax),rlonlo,rlonhi,rlatlo,rlathi,zlo,zhi,
     .       sqrt_nobs,rmjd0,rinterval,totalvariance,eofvariance,
     .       eofpercentage,cumper
      integer index(nmax),nvalid
      character*72 input,utgrid,ingrid,out,inmask
c
      nargs=iargc()
      call getarg(1,input)
      read (input,*) neof,rmjd0,rinterval
      call getarg(2,utgrid)
      lutgrid=lnblnk(utgrid)
      call getarg(3,inmask)
      maskx=0
      masky=mmax
      if (inmask.ne.'-') then
      call gridrd8(inmask,maskx,masky,mask,rlonlo,rlonhi,rlatlo,
     .         rlathi,zlo,zhi)
      if (maskx*masky.le.0)
     .   stop "error in input mask and/or mask size problems"
      else
      call getarg(4,inmask)
      call gridrd8(inmask,maskx,masky,mask,rlonlo,rlonhi,rlatlo,
     .         rlathi,zlo,zhi)
      do j=1,mmax
         if (mask(j).lt.1d20) mask(j)=0d0
      enddo
      endif
      nobs=nargs-3
      if (nobs.gt.nmax) stop "Too many input grids"
c
   13 format('reading ',a)
      nvalid=0
      totalvariance=0d0
      do k=1,nobs
         call getarg(k+3,ingrid)
         write (0,13) ingrid
         nx=maskx
         ny=masky
         call gridrd8(ingrid,nx,ny,grid,rlonlo,rlonhi,rlatlo,
     .            rlathi,zlo,zhi)
         if (nx*ny.ne.maskx*masky)
     .      stop "error in input grid and/or grid size problems"
         do i=1,nx*ny
            if (mask(i).lt.1d20) then
               nvalid=nvalid+1
	       if (nvalid.gt.amax) stop "Too many total observations"
               a(nvalid)=grid(i)-mask(i)
               totalvariance=totalvariance+(grid(i)-mask(i))**2
            endif
         enddo
      enddo
      totalvariance=totalvariance/nvalid
      ndim=nvalid/nobs
      call dsvdcmp(a,ndim,nobs,ndim,nobs,w,v)
c initialize index
      do k=1,nobs
         index(k)=k
      enddo
      call dqsort(index,w,nobs)
      sqrt_nobs=sqrt(dfloat(nobs))
      do i=1,neof
         write (out,101) utgrid(1:lutgrid),i
         k=index(nobs+1-i)
         nvalid=0
         do j=1,nx*ny
            if (mask(j).lt.1d20) then
               nvalid=nvalid+1
               grid(j)=a(nvalid+(k-1)*ndim)*w(k)/sqrt_nobs
            else
               grid(j)=1d35
            endif
         enddo
         write (6,*) 'EOF ',i,w(k),w(k)**2/ndim/nobs,totalvariance,
     .                      w(k)**2/ndim/nobs/totalvariance
         do j=1,nobs
	    write (6,*) rmjd0+(j-1)*rinterval,v(j+(k-1)*nobs)*sqrt_nobs
         enddo
      call gridwr8(out,nx,ny,grid,nx,rlonlo,rlonhi,rlatlo,
     .          rlathi)
      enddo
      cumper=0d0
      do i=1,nobs
         k=index(nobs+1-i)
         eofvariance=w(k)**2/ndim/nobs
         eofpercentage=eofvariance/totalvariance
         cumper=cumper+eofpercentage
         write (6,*) i,eofpercentage,cumper
      enddo
  101 format (a,'_eof',i2.2,'.grd')
      end
