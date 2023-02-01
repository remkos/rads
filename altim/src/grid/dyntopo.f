**DYNTOPO -- Compute X and Y velocities from dynamic topography
*+

*-
* 23-Feb-1998 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer mx,nx,ny,gridrd4,gridwr4,ios,iargc
      parameter (mx=801*801)
      real h(mx),vx(mx),vy(mx),x0,x1,y0,y1,z0,z1
      character*80 filenm

      if (iargc().ne.3) then
         write (0,600)
	 goto 9999
      endif

550   format (a)
600   format ('dyntopo -',
     |' Compute X and Y velocities from dynamic topography'//
     |'syntax:  dyntopo  dyntopo-grid velx-grid vely-grid')

* Read dyntopo grid

      call getarg(1,filenm)
      nx=0
      ny=mx
      ios=gridrd4(filenm,nx,ny,h,x0,x1,y0,y1,z0,z1)
      if (ios.ne.0) then
         write (0,550) 'dyntopo: error loading h grid'
	 goto 9999
      endif

* Compute Vx and Vy

      call dyntopo(h,vx,vy,nx,ny,nx,x0,x1,y0,y1)

* Write Vx and Vy grids

      call getarg(2,filenm)
      ios=gridwr4(filenm,nx,ny,vx,nx,x0,x1,y0,y1)
      if (ios.ne.0) then
         write (0,550) 'dyntopo: error writing vx grid'
	 goto 9999
      endif
      call getarg(3,filenm)
      ios=gridwr4(filenm,nx,ny,vy,nx,x0,x1,y0,y1)
      if (ios.ne.0) then
         write (0,550) 'dyntopo: error writing vy grid'
	 goto 9999
      endif

9999  end
