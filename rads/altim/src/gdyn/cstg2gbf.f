      program cstg2gbf

* This is an improved version of CSTG2GBF.
* Improved means: more accessible. If requires not to make any links
* to in and output files and thus is more user friendly.
*
* Syntax: cstg2gbf inputfile np.outputfile [ ql.outputfile ]
*
* If "ql.outputfile" it will be removed upon exit.
*
* This program requires files system.data and times.data to be found
* in directory $ALTIM/data/tables and cstg2gbf.nml in $ALTIM/nml
*
* 29-Jul-1998 - Remko Scharroo - Adopted from cstg2gbf by Ron Noomen
*  4-Aug-1998 - Splitted MAXPAS into MAXPNP and MAXPQL
*  6-Aug-1998 - Complete remodelling of the source
* 20-Jul-1999 - maxpass increased to 200000. Additional format checks.
* 10-Jan-2000 - removed check for year in century.
* 14-Jan-2000 - avoid possible rollover problem with epochs close to 24:00
* 17-Jan-2000 - no printing on checksum error. Create extra files with
*               updated data (on request: -u flag)
* 29-Jun-2000 - Convert '^M' and ' ' in CSTG input to '0'
* 30-Oct-2000 - Extended maxsat to 30
* 20-Feb-2004 - Added check on combination of k13 and k21 (Eelco)
*-----------------------------------------------------------------------
      include 'cstg2gbf.inc'

* GBF variables
      
      integer*2 j5(2)
      integer*4 idata(17),j1,j9,j13,j17,j37,j41,j53,refdate,iargc
      real*4    sigma,trop,iono,displ,cmas
      real*8	tfrac,range,date,date1/0d0/,deltat,co
      equivalence
     |  (idata(1),j1),    (idata(2),j5(1)),  (idata(3),j9),
     |  (idata(4),j13),   (idata(5),j17),    (idata(6),tfrac),
     |  (idata(8),range), (idata(10),j37),   (idata(11),j41),
     |  (idata(12),sigma),(idata(13),trop),  (idata(14),j53),
     |  (idata(15),iono), (idata(16),displ), (idata(17),cmas)

* CSTG variables

      integer*4 i1a,i1b,i13a,i13b,i25,i32,i37,i41,i48,iq25,iq30,iq34,
     |  iq37
      integer*4 k1,k8,k10,k13,k17,k19,k21,k25,k33,k39,k43,k44,k45,k46,
     |	k47,k48,k21hlp
      equivalence (i25,k48),(i32,iq25),(i37,iq30),(i41,iq34)

* Other variables

      integer*4 offset,itest,i,j,l,mdate,ih,it,ip,tdiff,statinfo,
     |	lnblnk
      character line*69,line0*69,dum*32
      logical   header,newpas,qlform/.false./,makeupdate/.false./

* Satellite variables

      integer*4 maxsat
      parameter (maxsat=30)
      integer*4 isat(maxsat)
      real*8	cmass(maxsat)
      namelist /cstg2gbf_nml/ isat,cmass

      integer*4 itwo(0:31),mjdy(00:99)
      integer*4	il(67),itot,isum,n1,n2,n3,nrec,nmeas(10,9999)

      real*8	vlight
      parameter (vlight=299792458.0d0)

* File handling

      integer*4	npout/0/,qlout/0/
      character*80 npfile
      logical   test

* Check arguments

      do i=1,iargc()
         call getarg(i,line)
         if (line(:2).eq.'-h') then
            write (0,600)
600         format (
     |'cstg2gbf - Convert SLR data in ONS format to GBF'//
     |'syntax: cstg2gbf [options] ',
     |'[np.outputfile [ql.outputfile]] <inputfile'//
     |'where [options] are:'/
     |'  -u : store the updates in separate files'/
     |'and:'/
     |'  inputfile    : input file of SLR data in ONS format'/
     |'    (input is read from standard input)'/
     |'  np.outputfile: output file of SLR normal points in GBF format'/
     |'  ql.outputfile: output file of SLR quick-look data in GBF ',
     |'format'/
     |'    (output files are optional: not generated when left blank)'//
     |'Further print output comes to standard output')
	    goto 9999
         else if (line(:2).eq.'-u') then
            makeupdate=.true.
         else if (npout.eq.0) then
            npout=20
            open (npout,file=line,form='unformatted')
	    npfile=line
         else
	    qlout=30
            open (qlout,file=line,form='unformatted')
         endif
      enddo

* Initialise variables

      npass=0
      nrec=0
      j5(1)=20
      j37=0
      j41=0
      iono=0
      displ=0
      header=.false.
      newpas=.false.
      n1=0
      n2=0
      n3=0
      i48=0

* Initialise arrays

      do j=1,9999
         do i=1,10
	    nmeas(i,j)=0
	 enddo
      enddo
      do i=1,maxpass
         p_t0(i)=1d30
	 p_t1(i)=-1d30
	 p_nr(i)=0
      enddo
      do i=00,99
         mjdy(i)=mdate(2,i*10000+0101)-1
      enddo	 
      itwo(0)=1
      do i=1,31
         itwo(i)=itwo(i-1)*2
      enddo

* Read namelist

      call checkenv('ALTIM',line,l)
      line(l+1:)='/nml/cstg2gbf.nml'
      open (9,file=line)
      read (9,nml=cstg2gbf_nml)
      
* CHECKING LINE-TYPES
*   Ni = 1 : Line is data-type identifier
*        2 : Line is a headerline
*        3 : Line is a dataline
*        4 : Line contains a wrong character
*        5 : Checksum of the line incorrect
* i=3 for the current line
* If HEADER=.FALSE. the next data-block will be ignored

* Reading lines from data-file

10    read (*,'(a69)',end=300) line
      if (line.eq.' ') goto 10
      nrec=nrec+1
      n1=n2
      n2=n3

* Detection Normal Point or Quicklook data

      if (line(1:5).eq.'88888') then
         qlform=.true.
         n3=1
         header=.true.
         goto 10
      else if (line(1:5).eq.'99999') then
         qlform=.false.
         n3=1
         header=.true.
         goto 10
      else if (.not.header) then
         goto 10
      endif

* Character check

      do i=1,69
	 j=ichar(line(i:i))
         if (line(i:i).eq.' ' .or. j.eq.13) then
	    line(i:i)='0'
         else if (line(i:i).ge.'0'.and.line(i:i).le.'9') then
	 else if (line(i:i).eq.'-') then
	 else
	    n3=4
	    if (n2.eq.1) then
	       write(*,'("Data-block will be ignored, there is a ",
     |"wrong character in the header-line:"/a69)') line
               header=.false.
	       goto 10
	    else if (n2.ge.2) then
	       write (*,'("Data-line will be ignored, there is a ",
     |"wrong character in the line:"/a69)') line   
	       goto 10
	    endif    
	 endif
      enddo	 

* Detection Headerline

      read (line(1:7),'(i7)') j
      do i=1,maxsat
         if (j.eq.isat(i)) goto 70
      enddo	 
      goto 200

* HEADER-LINE

70    n3=2
      header=.true.
      newpas=.true.
      line0=line
      do i=1,69
         if (line0(i:i).eq.'0') line0(i:i)=' '
      enddo

* Verification checksum header-line

      read (line0(1:54),'(bz,52i1,i2)') (il(i),i=1,52),isum
      itot=0
      do i=1,52
         itot=itot+il(i)
      enddo
      itot=mod(itot,100)
      if (itot.ne.isum.and.isum.ne.0) then
         write(*,'("the checksum is not correct in the header",
     |        "-line:",i3/a69/"process continues",
     |        " and checksum must be:",i3/)') isum,line,itot
      endif

* Reading data from headerline

      read (line0(1:51),100,err=101) k1,k8,k10,k13,k17,k19,
     |  k21,k25,k33,k39,k43,k44,k45,k46,k47,k48
100   format (bz,i7,i2,i3,i4,2i2,i4,i8,i6,i4,5i1,i4)
      goto 10
101   write (*,'("Error reading header line:"/a69)') line
      n3=4
      goto 10

* DATA-LINE

200   n3=3

* Checking if pass has a correct header-line

      if (n2.le.1) then
         write (*,'("Data-block ignored because of ",
     |         "missing header-line. Nrec:",i7)') nrec
         header=.false.
         goto 10
      else if (n2.ge.4 .and. n1.le.1) then
         write (*,'("Data-block ignored because of ",
     .         "error in header-line. Nrec:",i7)') nrec
         header=.false.
         goto 10
      else if (n1.ne.1 .and. n2.eq.2) then
         write (*,'("Data-block ignored because of ",
     .         "missing data-type identifier. Nrec:",i7)') nrec
         header=.false.
         goto 10
      else if (.not.header) then
         goto 10
      endif

* Mod by Eelco Doornbos
      IF(k13.eq.7810 .AND. k21.eq.8460) k13 = 6810
      IF(k13.eq.7405 .AND. k21.eq.8470) k13 = 6405

* Verification checksum data-line

      if (qlform) then
         read (line(1:69),'(67i1,i2)') (il(i),i=1,67),isum
         itot=0
         do i=1,67
            itot=itot+il(i)
         enddo
      else
         read (line(1:63),'(52i1,i2)') (il(i),i=1,52),isum
         itot=0
         do i=1,52
            itot=itot+il(i)
         enddo
      endif
      itot=mod(itot,100)
      if (itot.ne.isum .and. isum.ne.0) then
	 nmeas(10,k13)=nmeas(10,k13)+1
      endif

* Reading data from data-line

      if (qlform) then
         read (line(1:48),'(4i6,i5,i4,i3,i8)',err=290)
     |		i1a,i1b,i13a,i13b,iq25,iq30,iq34,iq37
         if (iq37.eq.0) iq37 = k25
      else
         read (line(1:48),'(4i6,i7,i5,i4,i3,4x,i1)',err=290)
     |		i1a,i1b,i13a,i13b,i25,i32,i37,i41,i48
      endif

* Calculate time in MJD since reference date

      date=dble(i1a)/86400d1+dble(i1b)/86400d7
      refdate=mjdy(k8)+k10
      
* Store date of first observation in block

      if (newpas) date1=date+refdate

* Check if 24 h. is crossed. If so, increase reference date by one.

      if (abs(date+refdate-date1).gt.0.5d0) refdate=refdate+1

* Time conversions

      if (k44.eq.7) then
         i=tdiff(date+refdate,deltat,co)
	 if (i.eq.1) then
	    nmeas(3,k13)=nmeas(3,k13)+1
	 else if (i.eq.3) then
	    nmeas(4,k13)=nmeas(4,k13)+1
	 endif   
         date=date-deltat/86400d0	! UTC(BIPM) -> UTC(USNO)
      else if (k44.eq.4) then
         i=tdiff(date+refdate,deltat,co)
	 if (i.eq.2) then
	    nmeas(5,k13)=nmeas(5,k13)+1
	 else if (i.eq.3) then
	    nmeas(6,k13)=nmeas(6,k13)+1
	 endif   
         date=date+(co-deltat)/86400d0	! UTC(GPS) -> UTC(USNO)
      else if (k44.eq.3) then
      					! No conversion
      else       					
        nmeas(2,k13)=nmeas(2,k13)+1	! Unknown
        goto 10
      endif

* The parameter J5(2) must be equal to 23 at this moment
*               J5(2)= 10*2 + K44 = 23        
* The factor 2 in '10*2' indicates the time of laser firing.
* (In MERIT-II format it used to be '10*1' to indicate the satellite time)
* At this stage the factor K44 must be equal to 3 to indicate the time has
* been transformed to UTC(USNO).                                        

      j5(2)=23
      if (i32.ne.0) then
        itest=5
      else

* No meteo data available; hence no tropospheric correction can be
* applied. Contrary to the description of the ONSITE data format,
* the observation is already corrected for this delay.

        itest=0
      endif

* Checking satellite

      offset=0
      cmas=0
      do i=1,maxsat
         if (k1.eq.isat(i)) then
	    cmas=cmass(i)
	    if (abs(cmas).lt.1e-10) offset=1
	    goto 220
	 endif 
      enddo	    
      nmeas(7,k13)=nmeas(7,k13)+1

220   j1=k1
      j9=k13
      j13=offset*itwo(28)+3*itwo(25)+itest*itwo(19)+1*itwo(18)
      j17=int(date+1d0)+(refdate-1)
      tfrac=(date+1d0)-int(date+1d0)
      range=dble(i13a)*1d6+dble(i13b)
      if (qlform) range=range-dble(iq37)
      range=range*0.5d-12*vlight+dble(cmas)
      sigma=dble(i25)*0.5d-12*vlight
      k21hlp=0
      if (statinfo(date+refdate,k13,dum,dum,dum,
     |    dum,k21hlp,dum,dum).gt.0) nmeas(8,k13)=nmeas(8,k13)+1

      if (k21.ne.0) then

* The next IF changes the units of the wavelength from 1.0 to 0.1 nm
* depending on the data value as follows:
* value 3000 to 9999: units 0.1 nm; value 0001 to 2999: units 1.0 nm. 

        if (k21.le.2999) k21=k21*10

        if ((nint(k21*0.1d0).ne.k21hlp).and.newpas) then
          write(*,296) k13,nint(dble(k21)*1.0d-1),k21hlp,nrec
 296      format('WARNING: Wavelength on data record for ',I4,
     |' <',I4,'> and system.data <',I4,'> mismatch (rec:',I7,')')
        endif
        k21hlp=nint(k21*0.1d0)
      endif

      trop=(k21hlp*1d-3+99d0)*1d-30
      ih=i41
      it=(i37+5)/10
      ip=(i32+5)/10
      j53=ih*itwo(24)+it*itwo(12)+ip

* Update data-array for making catalog
* Add 1 to number of passes for first measurement in a pass

      if (newpas) then
         npass=npass+1
	 if (npass.gt.maxpass) call fin('cstg2gbf: too many passes')
	 p_ql(npass)=(qlform .or. k43.eq.0)
         p_sat(npass)=k1
         p_sta(npass)=k13
         newpas =.false.
      endif
      p_t0 (npass)=min(p_t0(npass),date+refdate)
      p_t1 (npass)=max(p_t1(npass),date+refdate)
      p_nr (npass)=p_nr(npass)+1

* Write data to files

      if (p_ql(npass)) then
         if (qlout.gt.0) write (qlout) idata
      else if (npout.eq.0) then
      else if (makeupdate .and. i48.ne.0) then
         inquire (npout+i48,opened=test)
	 if (.not.test) then
	    l=lnblnk(npfile)
	    write (line,'(a,"_",i1)') npfile(:l),i48
	    open (npout+i48,file=line,form='unformatted')
	 endif
         write (npout+i48) idata
      else
         write (npout) idata
      endif	 	 
      goto 10

* Format error trap

290   write(*,'("Warning: the following line has a format error",
     |        "in characters 1-47:"/a69/"Line deleted"/)') line
      nmeas(9,k13)=nmeas(9,k13)+1
      goto 10

* Finishing up

300   write (*,'("Total number of lines read :",i7)') nrec
      do i=npout,qlout
         close (i)
      enddo
      call catalog
      call prerror(nmeas)
      write (*,'(/"Normal end of program cstg2gbf")')

9999  end
