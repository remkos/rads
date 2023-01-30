	program ascii2bimg
c-----------------------------------------------
c this program convert the ASCII distributed files for
c the Grenoble Tide Model to binary files in order to save
c disk space and time in reading the files.
c This Binary file format is interfaced with our graphical
c software.
c Jean-Marc Molines 15/11/1994
c------------------------------------------------
	parameter (nimax=720,njmax=361)
	real*4 wra(0:nimax,njmax), wrg(0:nimax,njmax)
	real a(30),G(30) 
	real*4 xmin,ymin,xma,yma,dx,dy,amask, Gmask
	integer ni,nj,narg
	character*80 filnama,filnamb,filheader
	character*80 line1, line2, line3, line4, line
	narg=iargc()
	iunit=10
c
	if (narg.ne.2) then
	  print *,
     |' >>>> ascii2bimg usage : ascii2bimg {file.model} {headerfile} '
	  print *,'     Binary file will be {file.model.bimg}'
	  print *,'   where ''model'' is the name of the model'
	  print *,'  headerfile is a file containing info on the model'
	  print *,'  as 3 records of 80 char. maximum.'
	  stop
	endif
c
	call getarg(1,filnama)
	call getarg(2,filheader)
	open(99,file=filheader,status='old',err=101)
	read(99,'(a)') line2
	read(99,'(a)') line3
	read(99,'(a)') line4
	close(99)
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
	      ii=k+(i-1)
	      if (ii.gt.720) print *,'ii > 720', ii
	      if (a(k).ne.amask.and.G(k).ne.Gmask) then
	        wra(ii,j)=a(k)
	        wrg(ii,j)=G(k)
	      else
	        wra(ii,j)=amask
	        wrg(ii,j)=amask
	      endif
	    enddo
	  enddo

c
c tide(-0.5,phi)=tide(359.5,phi)
c convenient for interpolation
c over a global solution
c
          wra(0,j)=wra(ni,j)
	  wrg(0,j)=wrg(ni,j)
	enddo	    
10	format(30f7.2)
20	format(30f7.1)

	close(iunit) 
 	l=lnblnk(filnama)
	filnamb=filnama(1:l)//'.bimg'
	   
	print *,' Writing  binary  file ..''will not be so long ...'
	open(11,file=filnamb,form='unformatted')
	z0=0.
	tim_tag=0.
	write(11)line1
	write(11)line2
	write(11)line3
	write(11)line4
	write(11)nimax+1,njmax,1,1,2,100
	write(11)xmin,ymin,dx,dy,amask
	write(11) z0
	write(11) tim_tag
	write(11) ((wra(i,j),i=0,nimax),j=1,njmax)
	write(11) ((wrg(i,j),i=0,nimax),j=1,njmax)
	close(11)
	goto 110
100	print *,' Problem opening ASCII file ... please check'
	goto 110
101	print *,' Problem opening header file ... please check'
110	continue
	end
	
	
 
