      program gridshift

      character arg*80,line*19
      integer	iarg,iargc
      real*8	shift/-0.145d0/,zmid

      if (iargc().eq.0) then
	 write (0,666)
	 goto 9999
666	 format (
     |  'shiftgrid: shift grid in vertical direction'//
     |  'syntax: gridshift [ options ] gridname(s)'//
     |  'where [ options ] are:'/
     |  ' shift=dist : add dist (cm) to heights (def: -14.5)')
      endif

      do iarg=1,iargc()
	 call getarg(iarg,arg)
	 if (arg(1:6).eq.'shift=') then
	    read (arg(7:),*) shift
	    shift=shift/100
	 else
      	    open (10,file=arg,status='old',form='unformatted',
     |	    	recl=19,access='direct')
      	    read (10,rec=1) line
      	    read (line,'(5x,e14.7)') zmid
      	    write (line,'("@GRID",e14.7)') zmid+shift
      	    write (10,rec=1) line
      	    close (10)
	 endif
      enddo

9999  end
