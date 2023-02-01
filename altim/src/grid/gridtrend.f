      implicit none
      integer nx,ny,mmax,nmax,nobs,nargs,iargc,k,l,mx,my
      parameter (mmax=361*181,nmax=100)
      real*4 grid(mmax*nmax),grida(mmax),gridb(mmax),gridr(mmax),
     .       rlonlo,rlonhi,rlatlo,rlathi,zlo,zhi
      real*8 x(nmax),y(nmax),a,b,r
      character*80 ingrid
c
      nargs=iargc()
      nobs=nargs
c
      mx=0
      my=mmax
      do k=1,nobs
         call getarg(k,ingrid)
         nx=mx
         ny=my
         call gridrd4(ingrid,nx,ny,grid(mx*my*(k-1)+1),
     |		rlonlo,rlonhi,rlatlo,rlathi,zlo,zhi)
         if (k.gt.1 .and. (nx.ne.mx .or. ny.ne.my)) stop "wrong size"
	 write (6,*) nx,ny,mx,my
	 mx=nx
	 my=ny
      enddo

      do l=1,mx*my
	 do k=1,nobs
	    x(k)=(k-1)*35e0/36525e0
	    y(k)=grid(mx*my*(k-1)+l)
	    if (y(k).gt.1e20) then
	       a=1e20
	       b=1e20
	       r=1e20
	       goto 20
	    endif
	 enddo
	 call regres(nobs,x,y,a,b,r)
20       continue
         grida(l)=a
         gridb(l)=b
         gridr(l)=r
      enddo

      call gridwr4('a.grd',nx,ny,grida,nx,
     |		rlonlo,rlonhi,rlatlo,rlathi)
      call gridwr4('b.grd',nx,ny,gridb,nx,
     |		rlonlo,rlonhi,rlatlo,rlathi)
      call gridwr4('r.grd',nx,ny,gridr,nx,
     |		rlonlo,rlonhi,rlatlo,rlathi)
      end
