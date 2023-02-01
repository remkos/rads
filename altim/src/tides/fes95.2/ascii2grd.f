	program ascii2grd
c-----------------------------------------------
c this program convert the ASCII distributed files for
c the Grenoble Tide Model to binary grid files in order to save
c disk space and time in reading the files.
c This Binary file format is interfaced with our graphical
c software.
c Remko Scharroo
c------------------------------------------------
	parameter (nimax=721,njmax=361)
	real*4 wra(nimax,njmax), wrg(nimax,njmax)
	real a(30),G(30) 
	real*4 xmin,ymin,xma,yma,dx,dy,amask, Gmask
	integer ni,nj,narg,iarg,iargc,l,m
	character*80 filnama,filnamb,fildir
	character*80 line1, line
	logical add,error
	real fact/1./

	narg=iargc()
	iunit=10
	rad=4*atan(1.)/180
c
	if (narg.le.0) then
	  print *,
     |' >>>> ascii2grd usage : ascii2grd {file.model(s)} '
	  print *,
     |		'Grid files will be {file.model.amp} and {file.model.pha}'
	  print *,
     |		'where ''model'' is the name of the model'
	  stop
	endif
c
	do iarg=1,narg
	call getarg(iarg,filnama)
	if (filnama(1:2).eq.'-a') then
	   add=.true.
	   write (6,*) 'adding FES95.2'
	   fact=10
	   write (6,*) 'multiplying phases by ',fact
	   goto 200
	else if (filnama(1:2).eq.'-e') then
	   error=.true.
	   write (6,*) 'generating error loading'
	   goto 200
	else if (filnama(1:2).eq.'f=') then
	   read (filnama(3:),*) fact
	   goto 200
	endif

	if (add) then
	   l=index(filnama,'.')
	   print *,'Loading old grid...'
	   nx=nimax
	   ny=njmax
	   fildir='/user/altim'
	   call checkenv('ALTIM',fildir,m)
	   filnamb=
     |fildir(:m)//'/data/grenoble/'//filnama(1:l)//'fes95.2.amp'
	   call gridrd4(filnamb,nx,ny,wra,xmin,xma,ymin,yma,zmin,zma)
	   filnamb=
     |fildir(:m)//'/data/grenoble/'//filnama(1:l)//'fes95.2.pha'
	   write (6,'(2i6,6f9.3)') nx,ny,xmin,xma,ymin,yma,zmin,zma
	   nx=nimax
	   ny=njmax
	   call gridrd4(filnamb,nx,ny,wrg,xmin,xma,ymin,yma,zmin,zma)
	   write (6,'(2i6,6f9.3)') nx,ny,xmin,xma,ymin,yma,zmin,zma
	else
	   do ny=1,njmax
	      do nx=1,nimax
		 wra(nx,ny)=0
		 wrg(nx,ny)=0
	      enddo
	   enddo
	endif

	open(iunit,file=filnama,status='old',err=100)
	print *,' Reading ASCII file ... be patient ...'

	read(iunit,'(a)') line
	read(line,*) xmin,xma
	line1=line(19:)
	print *,line1
	read(iunit,*) ymin,yma
	read(iunit,*) dx,dy
	read(iunit,*)   ni,nj
	read(iunit,*) amask,Gmask

	write(6,'(2f9.3)') xmin,xma
	write(6,'(2f9.3)') ymin,yma
	write(6,'(2f9.3)') dx,dy
	write(6,'(2i6)')   ni,nj
	write(6,'(2f9.3)') amask,Gmask

	do j=1,nj
	  do i=1,ni,30
	    read(iunit,10) (a(k),k=1,30)
	    read(iunit,10) (G(k),k=1,30)
	    do k=1,30
	      G(k)=G(k)*fact
	      ii=k+(i-1)
	      if (ii.gt.720) print *,'ii > 720', ii
	      if (a(k).ge.amask .or. G(k).ge.Gmask
     |	.or. wra(ii,j).gt.1e20) then
	        wra(ii,j)=1e30
	        wrg(ii,j)=1e30
	      else if (add) then
		x=wra(ii,j)*cos(wrg(ii,j)*rad)+a(k)*cos(G(k)*rad)
		y=wra(ii,j)*sin(wrg(ii,j)*rad)+a(k)*sin(G(k)*rad)
*		write (6,'(4f9.3,$)') wra(ii,j),wrg(ii,j),a(k),G(k)
	        wra(ii,j)=sqrt(x*x+y*y)
	        wrg(ii,j)=atan2(y,x)/rad
		if (wrg(ii,j).lt.0) wrg(ii,j)=wrg(ii,j)+360
*		write (6,'(2f9.3)') wra(ii,j),wrg(ii,j)
	      else if (error) then
		x=cos(G(k)*rad)-cos(G(k)*10*rad)
		y=sin(G(k)*rad)-sin(G(k)*10*rad)
	        wra(ii,j)=a(k)*sqrt(x*x+y*y)
	        wrg(ii,j)=atan2(y,x)/rad
	      else
		wra(ii,j)=a(k)
		wrg(ii,j)=G(k)
	      endif
	    enddo
	  enddo

          wra(ni+1,j)=wra(1,j)
	  wrg(ni+1,j)=wrg(1,j)
	enddo	    
10	format(30f7.2)
20	format(30f7.1)

	xma=xma+dx
	close(iunit) 
	print *,' Writing binary files ..''will not be so long ...'
	l=index(filnama,' ')-1
	filnamb=filnama(1:l)//'.amp'
	if (error) filnamb=filnama(1:l-4)//'error.amp'
	call gridwr4(filnamb,ni+1,nj,wra,nimax,xmin,xma,ymin,yma)
	filnamb=filnama(1:l)//'.pha'
	if (error) filnamb=filnama(1:l-4)//'error.pha'
	call gridwr4(filnamb,ni+1,nj,wrg,nimax,xmin,xma,ymin,yma)
200	enddo
	goto 110
100	print *,' Problem opening ASCII file ... please check'
	goto 110
101	print *,' Problem opening header file ... please check'
110	continue
	end
