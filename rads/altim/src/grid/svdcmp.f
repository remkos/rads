      program svdcmp

* This program makes a Singular Value Decomposition of several grids,
* equally seperated in time. Some people may call this an EOF analysis.
* The heart of the program is DSVDCMP from Numerical Recipes.
*
* Syntax:
*   svdcmp [ options ] grids ...
* (type just 'svdcmp' at the prompt for a full syntax)
*-
* 09-Oct-1998 -- Shaped up from SVD by Remko Scharroo
*-----------------------------------------------------------------------
      implicit none
      integer nx,ny,mmax,nmax,nobs,iargc,i,j,k,neof,amax,ndim
      parameter (mmax=361*181,nmax=400,amax=mmax*120)
      real*8 a(amax),v(nmax*nmax),w(nmax),grid(mmax),
     |       mask(mmax),t(nmax),x0,x1,y0,y1,z0,z1,
     |       sqrt_nobs,t0,dt,totalvariance,eofvariance,
     |       eofpercentage,cumper
      integer idx(nmax),nvalid,lnblnk
      character*80 arg,gridnm(nmax),pref,out,fmt

* Initialise

      nx=0
      ny=mmax
      neof=0
      nobs=0
      do j=1,mmax
         mask(j)=0d0
      enddo
      nvalid=0
      totalvariance=0d0
      pref='eof_'
      fmt='-'
      t0=0
      dt=1

* Scan arguments

      do i=1,iargc()
         call getarg(i,arg)
	 if (arg(:2).eq.'t=') then		! Set the start time and step
	    read (arg(3:),*,iostat=j) t0,dt
	 else if (arg(:4).eq.'fmt=') then
	    fmt=arg(5:)
	 else if (arg(:5).eq.'mask=') then	! Read a maskfile
            call gridrd8(arg(6:),nx,ny,mask,x0,x1,y0,y1,z0,z1)
      	    do j=1,nx*ny
               if (mask(j).lt.1d20) mask(j)=0d0
      	    enddo
	 else if (arg(:5).eq.'neof=') then
	    read (arg(6:),*) neof
	 else if (arg(:5).eq.'pref=') then
	    pref=arg(6:)
	 else
	    nobs=nobs+1
	    if (nobs.gt.nmax) call fin("Too many input grids")
	    gridnm(nobs)=arg
	 endif
      enddo

* Print Syntax

      if (nobs.eq.0) then
         write (0,600)
	 goto 9999
      endif
600   format ('svdcmp -- Do Singular Value Decomposition'//
     |'Syntax: svdcmp [ options ] grids ...'//
     |'where [options] are:'/
     |'  t=t0,dt   : Specify start time and time interval of grids',
     |' (def: t=0,1)'/
     |'  fmt=fmt   : Format for reading time out of file name',
     |' (example: ''fmt=(14x,f6.0)'')'/
     |'  neof=n    : Specify number of EOFs to produce',
     |' (def: equal to number of grids)'/
     |'  mask=file : Specify file to mask invalid points'/
     |'  pref=str  : Specify prefix for output files (def: ''eof_'')')

* If NEOF not specified, default it to NOBS

      if (neof.eq.0) neof=nobs

* If mask not specified, make them ourselves

      if (nx.eq.0) then
         do i=1,nobs
            call gridrd8(gridnm(i),nx,ny,grid,x0,x1,y0,y1,z0,z1)
	    if (nx*ny.eq.0) call fin("Error loading mask")
            do j=1,nx*ny
	       if (grid(j).gt.1d20) mask(j)=1d20
	    enddo
	 enddo
      endif

* Load all grids

      do i=1,nobs
         call gridrd8(gridnm(i),nx,ny,grid,x0,x1,y0,y1,z0,z1)
	 if (nx*ny.eq.0) call fin("Error loading grid")
	 if (fmt.eq.'-') then
	    t(i)=t0+(j-1)*dt
	 else
	    read (gridnm(i),fmt) t(i)
	 endif
         do j=1,nx*ny
	    if (mask(j).lt.1d20) then
               nvalid=nvalid+1
	       if (nvalid.gt.amax) call fin("Too many total observations")
               a(nvalid)=grid(j)
               totalvariance=totalvariance+grid(j)**2
            endif
         enddo
      enddo
      totalvariance=totalvariance/nvalid
      ndim=nvalid/nobs
      sqrt_nobs=sqrt(dfloat(nobs))

* Do the SVD

      call dsvdcmp(a,ndim,nobs,ndim,nobs,w,v)

* Initialise IDX, then do sort of indexes by value of W (power)

      do i=1,nobs
         idx(i)=i
      enddo
      call dqsort(idx,w,nobs)

* Write out the EOFs, starting with the highest power

      do i=1,neof
         k=idx(nobs+1-i)
         nvalid=0
         do j=1,nx*ny
            if (mask(j).lt.1d20) then
               nvalid=nvalid+1
               grid(j)=a(nvalid+(k-1)*ndim)*w(k)/sqrt_nobs
            else
               grid(j)=1d35
            endif
         enddo
         write (6,610) i,w(k),w(k)**2/ndim/nobs,totalvariance,
     |                      w(k)**2/ndim/nobs/totalvariance
         do j=1,nobs
	    write (6,620) t(j),v(j+(k-1)*nobs)*sqrt_nobs
         enddo
	 j=lnblnk(pref)
	 write (out,'(a,i3.3,".grd")') pref(:j),i
         call gridwr8(out,nx,ny,grid,nx,x0,x1,y0,y1)
      enddo
610   format ('EOF',i6,3d16.8,f12.8)
620   format (f15.3,f12.8)

* Write out EOF power percentages

      cumper=0d0
      do i=1,nobs
         k=idx(nobs+1-i)
         eofvariance=w(k)**2/ndim/nobs
         eofpercentage=eofvariance/totalvariance
         cumper=cumper+eofpercentage
         write (6,*) i,eofpercentage,cumper
      enddo
9999  end
