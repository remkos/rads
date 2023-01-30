      program xyz2deos

* Change XYZ data to DEOS grid file

      implicit none
      integer*4	nx,ny,kx,ky,i,maxgrd,ios,iargc
      real*8	x0/-180d0/,x1/180d0/,dx/1d0/,y0/-90d0/,y1/90d0/,dy/1d0/,
     |		x,y,z
      parameter (maxgrd=1000000)
      character*80 filenm/' '/,arg
      real*8	grid(maxgrd)

* Read arguments

      do i=1,iargc()
         call getarg(i,arg)
	 if (arg(:2).eq."x=") then
	    read (arg(3:),*,iostat=ios) x0,x1,dx
	 else if (arg(:2).eq."X=") then
	    read (arg(3:),*,iostat=ios) x0,x1,dx
	    x0=x0+dx/2
	    x1=x1-dx/2
	 else if (arg(:2).eq."y=") then
	    read (arg(3:),*,iostat=ios) y0,y1,dy
	 else if (arg(:2).eq."Y=") then
	    read (arg(3:),*,iostat=ios) y0,y1,dy
	    y0=y0+dy/2
	    y1=y1-dy/2
	 else
	    filenm=arg
	 endif
      enddo

* Check result

      if (filenm.eq.' ') then
         write (0,1300)
	 goto 9999
      endif
1300  format ('xyz2deos - Change XYZ data to DEOS grid file'//
     |'syntax: xyz2deos [options] gridfile < xyzfile'//
     |'where'//
     |'  xyzfile    : ASCII file with X,Y,Z triplets'/
     |'  gridfile   : DEOS grid file'//
     |'and [options] are:'//
     |'  x=x0,x1,dx :',
     |' set X-coordinate boundaries and resolution of nodes'/
     |'  X=x0,x1,dx :',
     |' set X-coordinate boundaries and resolution of cells'/
     |'  y=y0,y1,dy :',
     |' set Y-coordinate boundaries and resolution of nodes'/
     |'  Y=y0,y1,dy :',
     |' set Y-coordinate boundaries and resolution of cells')
1310  format ('xyz2deos: too many gridpoints ',i7,'>',i7)
      
      nx=nint((x1-x0)/dx+1)
      ny=nint((y1-y0)/dy+1)
      if (nx*ny.gt.maxgrd) then
         write (0,1310) nx*ny,maxgrd
	 goto 9999
      endif
      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)

* Initialise

      i=0
      do ky=1,ny
         do kx=1,nx
	    i=i+1
	    grid(i)=1d30
	 enddo
      enddo

* Read all the data and store values in matrices

10    read (*,*,end=100) x,y,z
      kx=nint((x-x0)/dx+1)
      ky=nint((y-y0)/dy+1)
      i=(ky-1)*nx+kx
      grid(i)=z
      goto 10

* Write out the grids

100   do i=1,4
         call gridwr8(filenm,nx,ny,grid,nx,x0,x1,y0,y1)
      enddo
9999  end
