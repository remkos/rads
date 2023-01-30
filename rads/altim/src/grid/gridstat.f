      implicit none

      integer*4	maxgrd
      parameter	(maxgrd=2000000)
      real*4	a(maxgrd),mean(maxgrd),rms(maxgrd),var(maxgrd)
      equivalence (a,var)

      integer*4	n(maxgrd),nmin/0/,ngrid/0/,nx,ny,mx,my,i,j,l,iargc,nval
      real*8	tmean,trms,tvar
      real*4	xmin,xmax,ymin,ymax,zmin,zmax

      character*80 arg,meanfile/' '/,rmsfile/' '/,varfile/' '/
      logical	options/.false./,files/.false./

      do i=1,iargc()
	 call getarg(i,arg)
	 if (arg(1:4).eq.'var=') then
	    varfile=arg(5:)
	    options=.true.
	 else if (arg(1:4).eq.'rms=') then
	    rmsfile=arg(5:)
	    options=.true.
	 else if (arg(1:5).eq.'mean=') then
	    meanfile=arg(6:)
	    options=.true.
	 else if (arg(1:4).eq.'min=') then
	    read (arg(5:),*) nmin
	 else
	    files=.true.
	 endif
      enddo

      if (.not.files .or. .not.options) then
	 write (0,600)
600	 format ('gridstat - compute statistical grid from various grids'//
     |  'syntax: gridstat [options] grid(s)'//
     |  'where [options] are:'/
     |  ' var=grid  : write variance to grid'/
     |  ' rms=grid  : write rms to grid'/
     |  ' mean=grid : write mean to grid'/
     |  ' min=n     : minimum nr to points to compute stats',
     |  ' (default = nr of grids)')
	 goto 9999
      endif

      mx=0
      my=maxgrd
      do j=1,maxgrd
         mean(j)=0
	 rms(j)=0
	 n(j)=0
      enddo

      do i=1,iargc()
	 call getarg(i,arg)
	 l=index(arg,' ')-1
	 if (index(arg,'=').eq.0) then
	    nx=mx
	    ny=my
	    ngrid=ngrid+1
	    call gridrd4(arg,nx,ny,a,xmin,xmax,ymin,ymax,zmin,zmax)
	    if (mx.ne.0 .and. (nx.ne.mx .or. ny.ne.my)) then
	       write (0,610) arg(:l)
610   format ('gridstat: file ',a,' not of same dimension as previous')
	    endif
	    mx=nx
	    my=ny
	    do j=1,nx*ny
	       if (a(j).lt.1d20) then
	          mean(j)=mean(j)+a(j)
	          rms(j)=rms(j)+a(j)**2
	          n(j)=n(j)+1
	       endif
	    enddo
	 endif
      enddo

      if (nmin.eq.0) nmin=ngrid

      do j=1,nx*ny
	 if (mean(j).gt.1d20 .or. n(j).lt.nmin) then
	    mean(j)=1d30
	    rms(j)=1d30
	    var(j)=1d30
	 else
	    mean(j)=mean(j)/n(j)
	    rms(j)=sqrt(rms(j)/n(j))
	    var(j)=sqrt(rms(j)**2-mean(j)**2)
	    tmean=tmean+mean(j)**2
	    trms=trms+rms(j)**2
	    tvar=tvar+var(j)**2
	    nval=nval+1
	 endif
      enddo

      tmean=sqrt(tmean/nval)
      trms=sqrt(trms/nval)
      tvar=sqrt(tvar/nval)
      write (6,700) ngrid,nmin,nx*ny,nval,tmean,trms,tvar
700   format (
     |'Number of grids         :',i9/
     |'Min nr of valids        :',i9/
     |'Total nr of grid points :',i9/
     |'Nr of valid grid points :',i9/
     |'Rms of mean grid        :',f9.4/
     |'Rms of rms grid         :',f9.4/
     |'Rms of variance grid    :',f9.4)


      if (varfile.ne.' ')
     |	call gridwr4(varfile,nx,ny,var,nx,xmin,xmax,ymin,ymax)
      if (rmsfile.ne.' ')
     |	call gridwr4(rmsfile,nx,ny,rms,nx,xmin,xmax,ymin,ymax)
      if (meanfile.ne.' ')
     |	call gridwr4(meanfile,nx,ny,mean,nx,xmin,xmax,ymin,ymax)

9999  end
