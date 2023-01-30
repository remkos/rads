      program gridamph

* Convert grids of C and S amplitude to Amplitude and Phase

      character arg*80
      integer	i,iargc,nx,ny,maxg,gridrd4
      parameter (maxg=721*181)
      real*4	cgrid(maxg),sgrid(maxg),agrid(maxg),pgrid(maxg),
     |		x0,x1,y0,y1,z0,z1,rad

      rad=atan(1.)/45.

      if (iargc().ne.4) then
	 write (0,666)
	 goto 9999
666	 format (
     |  'gridamph: convert grids of C and S amplitude to Amplitude and',
     |  ' Phase'//
     |  'syntax: gridamph C-grid S-grid Amp-grid Phase-grid')
      endif

      call getarg(1,arg)
      nx=0
      ny=maxg
      i=gridrd4(arg,nx,ny,cgrid,x0,x1,y0,y1,z0,z1)
      if (i.ne.0 .or. nx*ny.le.0) stop

      call getarg(2,arg)
      i=gridrd4(arg,nx,ny,sgrid,x0,x1,y0,y1,z0,z1)
      if (i.ne.0 .or. nx*ny.le.0) stop

      do i=1,nx*ny
         agrid(i)=sqrt(cgrid(i)**2+sgrid(i)**2)
	 pgrid(i)=atan2(sgrid(i),cgrid(i))/rad
	 if (pgrid(i).lt.0) pgrid(i)=pgrid(i)+360
      enddo

      call getarg(3,arg)
      call gridwr4(arg,nx,ny,agrid,nx,x0,x1,y0,y1)
      call getarg(4,arg)
      call gridwr4(arg,nx,ny,pgrid,nx,x0,x1,y0,y1)

9999  end
