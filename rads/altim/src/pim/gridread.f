**GRIDREAD -- Read grid or combination of grids
*+
      SUBROUTINE GRIDREAD (GRID, WK1, WK2, NX, NY, FILE,
     |  XMIN, XMAX, YMIN, YMAX)
      REAL GRID(*), WK1(*), WK2(*),
     |  XMIN, XMAX, YMIN, YMAX
      INTEGER NX, NY
      CHARACTER*(*) FILE
      include "pim.inc"
*
* NAME may be a construction of filenames separated by '-' or '+' to
* subtract or add grids
*
      character*1024 string,filenm
      integer i,l,nx0,ny0,nx1,ny1,n,npnt,mode
      real z,xmin1,xmax1,ymin1,ymax1,rad
      logical first

      first=.true.
      string=file
*     nx=0
*     ny=maxgrd
      nstack=1
      mode=11
      rad=atan(1.)/45

10    if (string.eq.' ') goto 100
      call strip(string,filenm)
      l=index(filenm,' ')-1

      if (filenm.eq.'+') then
	 mode=1
	 goto 10
      else if (filenm.eq.'-') then
	 mode=2
	 goto 10
      else if (filenm.eq.'*') then
	 mode=3
	 goto 10
      else if (filenm.eq.'amp{') then
	 mode=4
	 goto 10
      else if (filenm.eq.'pha{') then
	 mode=5
	 goto 10
      else if (filenm.eq.'re{') then
	 mode=6
	 goto 10
      else if (filenm.eq.'im{') then
	 mode=7
	 goto 10
      else if (filenm.eq.'&') then
	 mode=11
	 goto 10
      else if (filenm.eq.'ave{') then
	 mode=11
	 goto 10
      else if (filenm.eq.'rms{') then
	 mode=12
	 goto 10
      else if (filenm.eq.'std{') then
	 mode=13
	 goto 10
      else if (filenm.eq.'int{') then
	 mode=14
	 goto 10
      else if (filenm.eq.'slp{') then
	 mode=15
	 goto 10
      else if (filenm.eq.'cor{') then
	 mode=16
	 goto 10
      else if (filenm.eq.'fit{') then
	 mode=17
	 goto 10
      else if (filenm.eq.'sum{') then
	 mode=18
	 goto 10
      else if (filenm.eq.'}') then
         if (mode.gt.10) then
            call copy(npnt,grid,wk1)
            call gridstck(npnt,nstack,tx,wk1,wk2)
            call copy(npnt,wk2((mode-11)*npnt+1),grid)
	 endif
	 goto 10
      endif

      if (first) then
	 call gridread2(filenm(:l),nx,ny,wk1,xmin,xmax,ymin,ymax)
*	 if (nx.ne.nx0 .or. ny.ne.ny0) stop 'something is wrong'
         npnt=nx*ny
	 tx(1)=35.0/2
      else
	 xmin1=xmin
	 xmax1=xmax
	 ymin1=ymin
	 ymax1=ymax
	 nx1=0
	 ny1=maxgrd
	 call gridread2(filenm(:l),nx1,ny1,wk1,xmin1,xmax1,ymin1,ymax1)
	 if (abs(xmin1-xmin).gt.1e3.or.abs(xmax1-xmax).gt.1e3
     |		.or.abs(ymin1-ymin).gt.1e3.or.abs(ymax1-ymax).gt.1e3) then
	     write (6,*) xmin,xmax,ymin,ymax
	     write (6,*) xmin1,xmax1,ymin1,ymax1
     		stop "gridread: trying to combine incompatible grids"
	 endif
	 if (nx1.ne.nx .or. ny1.ne.ny) then
	    write (0,550) 'Adjusting resolution of '//filenm(:l)//' ...'
	    call grdinter(wk1,nx1,ny1,wk2,nx,ny)
	    call copy(npnt,wk2,wk1)
	 endif
      endif
      if (first) then
	 call copy(npnt,wk1,grid)
	 first=.false.
      else if (mode.eq.1) then
	 write (0,550) 'Adding grids ...'
	 do i=1,npnt
	    grid(i)=grid(i)+wk1(i)
	    if (grid(i).gt.1e20 .or. wk1(i).gt.1e20) grid(i)=1e35
	 enddo
      else if (mode.eq.2) then
	 write (0,550) 'Subtracting grids ...'
	 do i=1,npnt
	    grid(i)=grid(i)-wk1(i)
	    if (grid(i).gt.1e20 .or. wk1(i).gt.1e20) grid(i)=1e35
	 enddo
      else if (mode.eq.3) then
	 write (0,550) 'Multiplying grids ...'
	 do i=1,npnt
	    grid(i)=grid(i)*wk1(i)
	    if (grid(i).gt.1e20 .or. wk1(i).gt.1e20) grid(i)=1e35
	 enddo
      else if (mode.eq.4) then
	 write (0,550) 'Compute amplitude of complex grids ...'
	 do i=1,npnt
	    grid(i)=sqrt(grid(i)**2+wk1(i)**2)
	    if (grid(i).gt.1e20 .or. wk1(i).gt.1e20) grid(i)=1e35
	 enddo
      else if (mode.eq.5) then
	 write (0,550) 'Compute phase of complex grids ...'
	 do i=1,npnt
	    grid(i)=atan2(wk1(i),grid(i))/rad
	    if (grid(i).gt.1e20 .or. wk1(i).gt.1e20) grid(i)=1e35
	 enddo
      else if (mode.eq.6) then
	 write (0,550) 'Compute real component of amp/pha grids ...'
	 do i=1,npnt
	    grid(i)=grid(i)*cos(wk1(i)*rad)
	    if (grid(i).gt.1e20 .or. wk1(i).gt.1e20) grid(i)=1e35
	 enddo
      else if (mode.eq.7) then
	 write (0,550) 'Compute imaginary component of amp/pha grids ...'
	 do i=1,npnt
	    grid(i)=grid(i)*sin(wk1(i)*rad)
	    if (grid(i).gt.1e20 .or. wk1(i).gt.1e20) grid(i)=1e35
	 enddo
      else if (mode.gt.10) then
	 write (0,550) 'Stacking grids ...'
	 if ((nstack+1)*npnt.gt.maxgrd) stop "grid too large"
	 call copy(npnt,wk1,wk1(nstack*npnt+1))
	 nstack=nstack+1
	 tx(nstack)=tx(nstack-1)+2*tx(1)
      endif
      goto 10

100   continue

* Eliminate bad points

      do i=1,npnt
         if (grid(i).lt.rminb .or. grid(i).gt.rmaxb) grid(i)=1e35
      enddo

550   format (a)
      end

      subroutine gridread2(filenm,nx,ny,work,xmin,xmax,ymin,ymax)
      character*(*) filenm
      integer nx,ny,l,i,sgridrd4
      real work(*),xmin,xmax,ymin,ymax,zmin,zmax,z

      l=len(filenm)
      if (filenm(1:1).eq.'(' .and. filenm(l:l).eq.')') then
	 write (0,550) 'Generating grid ',filenm,' ...'
	 read (filenm(2:l-1),*) z
	 nx=2
	 ny=2
	 do i=1,nx*ny
	    work(i)=z
	 enddo
	 zmin=z
	 zmax=z
      else
 	 write (0,550) 'Reading ',filenm,' ...'
         i=sgridrd4(filenm,nx,ny,work,xmin,xmax,ymin,ymax,zmin,zmax)
         if (i.ne.0) stop "gridread: error loading grid"
      endif
550   format(a,a,a)
      end
