**LGRID - list contents of binary grid in ASCII format
*+
      program lgrid
*
* Program to list the DEOS binary grid files
*
* syntax:  lgrid [options] gridname
* 
* with optional arguments:
*   fmt=fmt  : specify format for values, default=10F8.3 
*   inv=inv  : specify invalid value, default= 999.999  
*   -xyz     : write as XYZ (default stream output)
*
* required arguments:
*   gridname : filename of the grid
* 
* Output format (to standard output):
* Header: (6F10.4,2I6) XMIN,XMAX,YMIN,YMAX,DX,DY,NX,NY
* Body  : ((10F8.3 )) (GRID(K),K=1,NX*NY)
*-
* 26-Nov-1997 - Distributed version - Remko Scharroo (DUT/DEOS)
*  8-Nov-2003 - Added XYZ output
*-----------------------------------------------------------------------

      implicit none
      integer	mx,nx,ny,k,iargc,kx,ky
      parameter (mx=4000000)
      real*4    grid(mx),xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,
     |		inv/999.999/
      logical	xyz/.false./
      character*80 filenm/' '/,arg,fmt*7/'10F8.3 '/,inval*10

* Scan arguments
      
      do k=1,iargc()
         call getarg(k,arg)
	 if (arg(:4).eq.'fmt=') then
	    fmt=arg(5:)
	 else if (arg(:4).eq.'inv=') then
	    read (arg(5:),*) inv
	 else if (arg(:4).eq.'-xyz') then
	    xyz=.true.
	 else
	    filenm=arg
	 endif
      enddo

* Test arguments

      if (filenm.eq.' ') then
	 write (inval,'('//fmt//')') inv
         write (0,600) fmt,inval,fmt
	 goto 9999
      endif

* Load grid

      nx=0
      ny=mx
      call gridrd4(filenm,nx,ny,grid,xmin,xmax,ymin,ymax,zmin,zmax)
      dx=(xmax-xmin)/(nx-1)
      dy=(ymax-ymin)/(ny-1)

* XYZ output

      if (xyz) then
	 k=0
         do ky=1,ny
            do kx=1,nx
	       k=k+1
	       if (grid(k).lt.1e20)
     |		write (*,'('//fmt//')') xmin+(kx-1)*dx,ymin+(ky-1)*dy,grid(k)
	    enddo
         enddo

* Data stream output: first convert the invalid values

      else

         do k=1,nx*ny
	    if (grid(k).gt.1e20) grid(k)=inv
         enddo
         write (*,'(6f10.4,2i6)') xmin,xmax,ymin,ymax,dx*60,dy*60,nx,ny
         write (*,'('//fmt//')') (grid(k),k=1,nx*ny)
      endif

* Syntax

600   format('lgrid: list contents of grid in ASCII format'//
     |'syntax:  lgrid [options] gridname'//
     |'with optional arguments:'/
     |'  fmt=fmt  : specify format for values, default=',a/
     |'  inv=inv  : specify invalid value, default=',a/
     |'  -xyz     : print as XYZ (default stream output)'//
     |'required arguments:'/
     |'  gridname : filename of the grid'//
     |'Stream output format (to standard output):'/
     |'Header: (6F10.4,2I6) XMIN,XMAX,YMIN,YMAX,DX,DY,NX,NY'/
     |'Body  : ((',a,')) (GRID(K),K=1,NX*NY)'//
     |'XYZ output format (to standard output):'/
     |'X Y Z (starting in lower left corner)')
9999  end
