      program grid2xtra

* Convert DUT/SSR&T GRID format to KNMI EXTRA format

      implicit none
      integer nargs,nobs,iargc,mx,my,nx,ny,k,mt,ivar,layer,length,
     .        i,j,gridrd8
      parameter (mx=180,my=60,mt=160)
      real*8 rlonlo,rlonhi,rlatlo,rlathi,zlo,zhi
      character*72 ingrid
      real*8 grid(mx*my*mt)
      nargs=iargc()
      nobs=nargs-1
      call getarg(nargs,ingrid)
      open (10,file=ingrid,status='new',form='unformatted')
      j=1
      do k=1,nobs
         call getarg(k,ingrid)
         nx=0
         ny=mx*my
         i=gridrd8(ingrid,nx,ny,grid(j),rlonlo,rlonhi,rlatlo,
     .            rlathi,zlo,zhi)
         if (nx*ny.le.0)
     .      stop "error in input grid and/or grid size problems"
         j=j+nx*ny
      enddo
      length=nx*ny
      j=1
      do i=1,nobs*length
         if (grid(i).gt.1d20) grid(i)=9d99
      enddo
      do k=1,nobs
         write (10) k,ivar,layer,length
         write (10) (grid(i),i=j,j+length-1)
         j=j+length
      enddo
      end
