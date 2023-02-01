**GDR2MAR -- Program to convert Envisat GDR to MAR files
*+
      program gdr2mar
*
* This program converts Envisat (FD/I)GDR to MAR files. It simply copies
* the relevant data from the lengthier GDR files to shorter MAR files.
* No conversions are made, no editing is applied, no bugs are fixed.
*
* The only things that are non-standard in the output MAR:
* - The output file name starts with "RA2_MAR" instead of "RA2_WWV"
* - The total size in the header field "TOT_SIZE" is not correct
* - The 1-Hz Ku-band and S-band peakiness (type: US, unit: 1d-3) are
*   copied into spare bytes 141-144 of the MAR
* - Only the altimeter data block is copied. Possible other data blocks
*   are stripped off.
*
* Syntax:
* gdr2mar <filenames>
*
* All files specified on the command line are converted. Original
* file names should be used. The output files will have the same name
* with the string _GDR_, _IGD_, or _FGD_ replaced by _MAR_
*-
* 24-Feb-2003 - 0302.0 - Created by Remko Scharroo
* 25-Feb-2003 - 0302.1 - Proper initialisation of mds_off
*  9-Mar-2003 - 0303.0 - Added syntax, verbose, srcdir=, destdir=
*-----------------------------------------------------------------------
      implicit none
      character*256 arg,srcdir/' '/,destdir/' '/
      integer*4	iarg,iargc,l,gdrfd,marfd,rflag,wflag,lnblnk,gdropen
      logical	verbose/.false./

* Print syntax when no argument is supplied

      if (iargc().eq.0) then
         write (*,1300)
	 goto 9999
      endif

* Initialise the use of FASTIO

      call ioconst(rflag,wflag,wflag)

* Scan the argument list for GDR files.
* Return error messages when the input file is not a GDR file or
* when the input or output file can not be opened.

      do iarg=1,iargc()
         call getarg(iarg,arg)
	 l=lnblnk(arg)
	 if (arg(:2).eq.'-v') then
	    verbose=.true.
	 else if (arg(:7).eq.'srcdir=') then
	    srcdir=arg(8:)
	 else if (arg(:8).eq.'destdir=') then
	    destdir=arg(9:)
	 else if (arg(l-2:l).ne.'.N1' .or. arg(l-61:l-58).ne.'RA2_')
     |		then
	    write (*,600) arg(:l)
	 else if (arg(l-57:l-54).eq.'GDR_' .or.
     |		  arg(l-57:l-54).eq.'IGD_' .or.
     |		  arg(l-57:l-54).eq.'FGD_') then
	    gdrfd=gdropen(rflag,'<-',srcdir,arg)
	    if (gdrfd.le.0) goto 100
	    arg(l-57:l-55)='MAR'
	    marfd=gdropen(wflag,'->',destdir,arg)
	    if (marfd.le.0) goto 100
	    call gdrconv(gdrfd,marfd)
	    call closef(gdrfd)
	    call closef(marfd)
	 else
	    write (*,600) arg(:l)
	 endif
100   enddo

* Formats

600   format (a,' is not an Envisat GDR file')
1300  format ('gdr2mar - Convert Envisat GDR to MAR files'//
     |'syntax: gdr2mar [options] gdr-filename ...'//
     |'where [options] are:'/
     |'  -v          : be verbose'/
     |'  srcdir=dir  : specify source (GDR) directory (default: .)'/
     |'  destdir=dir : specify destination (MAR) directory',
     |' (default: .)')
9999  end

************************************************************************

      function gdropen(flag,arrow,dir,file)

* Open GDR file

      integer*4 gdropen,flag,l,openf,lnblnk
      character*(*) dir,file,arrow
      character*256 path

* Extend the file name to path name (include directory name)

      if (dir.eq.' ') then
         path=file
      else
         l=lnblnk(dir)
	 path=dir
	 path(l+1:l+1)='/'
	 path(l+2:)=file
      endif
      l=lnblnk(path)

* Open file for use with FASTIO routines

      write (*,600) arrow,path(:l)
      gdropen=openf(flag,o'666',path)
      if (gdropen.le.0) write (*,610) path(:l)

* Formats

600   format (a,' ',a)
610   format ('   Error opening file ',a)
      end

************************************************************************

      subroutine gdrconv(fdin,fdout)

* Convert the GDR file attached to fdin to a MAR file at fdout

      implicit none
      integer*4 fdin,fdout,headcnt,i,j,k,l,m,readf,writef
      integer*4 mds_off,mds_dim,mds_num_rec,mds_rec_dim,lnblnk
      character line*512,id*40,value*80
      logical	ra2dsd
      byte	app(2492)

* Scan the headers, take out the important information and change them

      headcnt=0
      ra2dsd=.false.
      mds_rec_dim=2492
      mds_num_rec=0
      mds_off=999999999

10    m=readf(fdin,0,line)
      headcnt=headcnt+m
      l=m-1
      if (line(:l).eq.' ') goto 100
      i=index(line(:l),'"')
      j=index(line(:l),'<')
      k=index(line(:l),'=')
      l=lnblnk(line(:l))
      if (i.gt.0) then
        value=line(i+1:l-1)	! string
      else if (j.gt.0) then
         value=line(k+1:j-1)	! value with unit
      else
         value=line(k+1:l)	! value without unit
      endif
      id=line(:k-1)

      if (id.eq.'PRODUCT') then
         line(14:16)='MAR'
      else if (id.eq.'SPH_DESCRIPTOR') then
         line(25:28)='IMAR'
      else if (id.eq.'DS_NAME') then
         if (value.eq.'RA2_DATA_SET_FOR_LEVEL_2') then
	    line(10:37)='RA2_OCEAN_DATA_FOR_LEVEL_2'
	    ra2dsd=.true.
	 else
	    ra2dsd=.false.
	 endif
      else if (id.eq.'DS_OFFSET') then
         if (ra2dsd) then
	    read (value,*) mds_off
	 else
	    write (line(12:31),600) 0
	 endif
      else if (id.eq.'DS_SIZE') then
         if (ra2dsd) then
	    read (value,*) mds_dim
	    mds_num_rec=mds_dim/mds_rec_dim
	    write (line(10:29),600) mds_num_rec*356
	 else
	    write (line(10:29),600) 0
	 endif
      else if (id.eq.'NUM_DSR') then
         if (.not.ra2dsd) write (line(10:19),610) 0
      else if (id.eq.'DSR_SIZE') then
         if (ra2dsd) then
	    write (line(11:20),610) 356
         else
	    write (line(11:20),610) 0
	 endif
      endif
100   k=writef(fdout,m,line)
      if (headcnt.lt.mds_off) goto 10

* All headers are processed. Start reading the data block.
* Now copy all the relevant information out of the GDR file into
* the MAR file. Also copy the Ku- and S-band peakiness into the spare
* bytes 141-144 of the MAR file.

      do i=1,mds_num_rec
         l=readf(fdin,mds_rec_dim,app)
	 k=writef(fdout, 40,app(   1))		!   1- 40
	 k=writef(fdout,  8,app(  81))		!  41- 48
	 k=writef(fdout,  8,app( 301))		!  49- 56
	 k=writef(fdout, 16,app( 469))		!  57- 72
	 k=writef(fdout, 56,app(1205))		!  73-128
	 k=writef(fdout, 12,app(1597))		! 129-140
*	 k=writef(fdout, 56,app(1889))		! 141-196
	 k=writef(fdout,  4,app(2473))		! 141-144 (peakiness)
	 k=writef(fdout, 52,app(1893))		! 145-196
	 k=writef(fdout,148,app(2305))		! 197-344
	 k=writef(fdout, 12,app(2477))		! 345-356
      enddo
600   format (i20.20)
610   format (i10.10)
      end
