      PROGRAM lres

* This program computes and prints the statistics of satellite tracking
* data stored in GEODYN II residual file format. It can also handle the
* .res files created by lresconv.
*
* Lres also has in editing mode that allows to screen SLR data using one
* of the interactive devices supported by PGPLOT (usually this is an X-windows
* terminal).
*
* For details of the use of this program, read the file "lres.man".
*
* Supported tracking data types are:
* - SLR
* - Altimetric heights
* - Altimeter crossovers
* - DORIS
* - PRARE range and range rate
*
* Latest changes:
* 29-Dec-1996 - Print extra mark for high elevation passes.
* 24-Aug-1999 - Minor upgrades
* 26-Oct-1999 - Add fitedit option and resfitedit routine
*  2-Nov-2000 - Add -D option to show statistics of deleted measurements
*  8-Nov-2000 - Rearrange effect of -D option
* 14-Nov-2000 - Sort PRARE data properly.
*               Measurement type added to lres.del
* 08-Mar-2002 - ED: No interactive editing for DOPPLER data
* 08-Jul-2002 - ED: Print satids with leading zero (e.g. 0200901 for Envisat)
* 18-Oct-2002 - ED: Initialize ntypes and nsats and use iargc() for compatibility with Intel f95 compiler
*-
* Created and maintained by Remko Scharroo
*-----------------------------------------------------------------------
      implicit none 

* Declaration of variables

      integer maxobs,maxpas
      parameter(maxobs=400,maxpas=5000)
      real*8 t1985
      parameter(t1985=(46066-30000)*86400d0)
      real*8 times(maxobs,maxpas),resid(maxobs,maxpas),
     |sigma(maxobs,maxpas),deriv(maxobs,maxpas),elev(maxobs,maxpas)
      integer ipass(4,maxpas),isort(maxpas,2),noise(maxpas),
     |itype,nsta,pobs(maxpas),kpass,ie0,ie1
      real*8 tmin,tmax,emax,
     |prms(maxpas),pmean(maxpas),pelev(maxpas)
      character cpass(maxpas)*8,prognm*80,
     |date0*15,date1*15,filenm*80,typenm1(100)*8,typenm2(100)*32,
     |mark*2,sitename(maxpas)*32,delfile*80
      integer maxtypes
      parameter (maxtypes=10)
      integer itypes(maxtypes),ntypes,isats(maxtypes),nsats
      logical timesort,zero,altim,long,exist,fitedit,delstat
      integer iarg,nfile,npass,ip,iplast/0/,iobs,ista,isat,nobs,i,j,l,
     |pgbeg,idel0,idel1,iter,iobs_old,maxiter,highelev(2), iargc
      character dev*80,line(maxpas)*120,text*120
      integer lnblnk
      real*8 rmean,rrms,a,b,fit,dum(3),editmult
      character*80 arg

* Common blocks

      include "lres.inc"

* Namelist

      namelist /lres_nml/ dev,rgb,xmin,xmax,ymin,ymax,long,altim,
     |	fitedit,timesort,maxiter,editmult,highelev

* Initialisations

      data delfile/' '/,zero/.true./,dev/' '/,delstat/.false./
      ntypes = 0 
      nsats = 0
      nsta = 0

* Open namelist and get defaults. Check for files called 'lres.nml'
* in the directory $ALTIM/nml, $HOME/nml, and the current directory

      text='/user/altim'
      call checkenv('ALTIM',text,l)
      open(7,file=text(:l)//'/nml/lres.nml')
      read(7,lres_nml)
      close(7)
      call checkenv('HOME',text,l)
      inquire(file=text(:l)//'/nml/lres.nml',exist=exist)
      if (exist) then
         open(7,file=text(:l)//'/nml/lres.nml')
         read(7,lres_nml)
         close(7)
      endif
      inquire(file='lres.nml',exist=exist)
      if (exist) then
         open(7,file='lres.nml')
         read(7,lres_nml)
         close(7)
      endif

* More initialisations

      call getarg(0,prognm)
      iarg=0
      filenm='fort.19'
      nfile=0
      typenm1( 51)='     SLR'
      typenm1( 40)='   DORIS'
      typenm1( 45)='PRARERAN'
      typenm1( 46)='PRAREDOP'
      typenm1( 99)='  HEIGHT'
      typenm1(100)='   XOVER'
      typenm2( 51)='SLR range'
      typenm2( 40)='DORIS doppler'
      typenm2( 45)='PRARE range'
      typenm2( 46)='PRARE doppler'
      typenm2( 99)='Altimeter height'
      typenm2(100)='Altimeter xovers'

* Scan the arguments

      do iarg=1,iargc()
      call getarg(iarg,arg)
      if (arg.eq.'-t') then
	 timesort=.true.
      else if (arg.eq.'-s') then
	 timesort=.false.
      else if (arg.eq.'-z') then
	 zero=.false.
      else if (arg.eq.'-l') then
	 long=.true.
      else if (arg.eq.'-a') then
	 altim=.false.
      else if (arg(1:2).eq.'-d') then
	 if (arg(3:).eq.' ') then
	    delfile='lres.del'
	 else
	    delfile=arg(3:)
	 endif
      else if (arg(1:2).eq.'-e') then
	 dev=arg(3:)
	 if (dev.eq.' ') dev='/xs'
	 if (delfile.eq.' ') delfile='lres.del'
	 altim=.false.
	 long=.true.
	 timesort=.true.
      else if (arg(1:2).eq.'-f') then
         fitedit=.true.
	 altim=.false.
      else if (arg(1:2).eq.'-D') then
	 delstat=.true.
      else if (arg(1:1).eq.'-') then
	 write (*,600)
	 goto 9999
      else
         filenm = arg
      endif
      enddo

* Open residual file and read it

      nfile=nfile+1
      call resread(filenm,maxobs,maxpas,ipass,cpass,
     |	times,resid,sigma,deriv,elev)

* Sort the passes in order of the station ID and satellite ID
* Use a little trick to get the PRARE data sorted correctly

      npass=0
      do ip=1,maxpas
	 if (ipass(3,ip).gt.0) then
	    npass=npass+1
            isort(npass,1)=ip
	    isort(npass,2)=
     |		ipass(1,ip)*100000+ipass(2,ip)/100
	    if (ipass(4,ip).eq.46) isort(npass,2)=isort(npass,2)+1
         endif
      enddo
      call bubble(isort(1,1),isort(1,2),npass)

* Process each pass one by one

      select=.false.
      do ii=1,npass
	 ip=ii
	 if (.not.timesort) ip=isort(ii,1)

* Copy pass buffer values to local variables

	 ista=ipass(1,ip)
	 isat=ipass(2,ip)
         nobs=ipass(3,ip)
	 itype=ipass(4,ip)
	 mark=' '

* Store data type and satellite number in respective buffers

	 call types (itypes,maxtypes,ntypes,itype)
	 call types (isats,maxtypes,nsats,isat)

* Determine start and end of pass

	 tmin=1d20
	 tmax=-1d20
         do i=1,nobs
	    tmin=min(tmin,times(i,ip))
	    tmax=max(tmax,times(i,ip))
	 enddo
	 if (nobs.gt.0) then
            call strf1985(date0,'%y%m%d %H:%M:%S',nint(tmin-t1985))
            call strf1985(date1,'%y%m%d %H:%M:%S',nint(tmax-t1985))
	 else
	    date0=' '
	    date1=' '
	 endif

* Eliminate low elevation points

	 do i=1,nobs
	    if (elev(i,ip).lt.10d0) sigma(i,ip)=0
	 enddo

* If statistics of the deleted measurements are requested, change the
* sign of the sigma

	 if (delstat) then
	    do i=1,nobs
	       sigma(i,ip)=-sigma(i,ip)
	       if (sigma(i,ip).eq.0d0) sigma(i,ip)=1
	    enddo
	 endif

* We have a problem with some PRARE-DOPPLER data.
* Attempt to overcome:
*
*        do i=1,nobs
*	    if (itype.eq.46 .and. abs(resid(i,ip)).gt.3d-3)
*    |		sigma(i,ip)=-abs(sigma(i,ip))
*        enddo

* Compute the statistics of the observation residuals for this pass.
* Computed are: mean, rms, fit (linear regression in range-rangerate
* diagram), and rms of fit.
* Only valid measurements and those above 10 degrees elevation are
* used.

	 iter=1
	 iobs_old=0
100	 call resfit(nobs,deriv(1,ip),resid(1,ip),sigma(1,ip),
     |		iobs,rmean,rrms,dum(1),a,b,fit)

* If this option is selected, do an editing based on the data fit
* Only do this when there are more then two points

         if (fitedit .and. iobs.ne.iobs_old .and. iobs.gt.2
     |		.and. iter.le.maxiter) then
	    iobs_old=iobs
	    call resfitedit(nobs,deriv(1,ip),resid(1,ip),sigma(1,ip),
     |		a,b,fit,editmult)
	    goto 100
	 endif

* Determine if pass if of high elevation

	 emax=0d0
	 ie0=0
	 ie1=0
	 do i=1,nobs
	    if (sigma(i,ip).gt.0) then
	       if (elev(i,ip).gt.emax) then
	          ie0=ie0+1
		  emax=elev(i,ip)
	       else
	          ie1=ie1+1
	       endif
	    endif
	 enddo
         if (emax.gt.highelev(1) .and. min(ie0,ie1).ge.highelev(2)) then
	    mark(1:1)='*'
	 endif

* Store pass residual statistics, pass dates and site name

	 pmean(ip)=rmean
	 prms(ip)=rrms
	 pobs(ip)=iobs
	 pelev(ip)=emax
	 if (ipass(4,ip).eq.99 .or. ipass(4,ip).eq.100) then
	    sitename(ip)='Global oceans'
	 else if (nobs.gt.0) then
	    call statinfo(tmin/86400d0+30000d0,ista,
     |		sitename(ip),dum,dum,noise(ip),dum,dum,dum)
	 endif

* Determine if pass is conspicuous

	 if ((iobs.le.2 .and. rrms.lt.1d-4) .or. iobs.ne.nobs .or.
     |		fit*1d2.gt.max(noise(ip),4) .or.
     |		abs(a/2).gt.sigma(1,ip) .or. abs(b).gt.50e-6) then
*	    write(*,*) iobs,rrms,nobs,fit,noise(ip),a,sigma(1,ip),b
	    mark(2:2)='!'
	 endif

* Print statistics per pass

	 line(ip)=' '
	 if (zero.or.iobs.gt.0) then
           if (itype.eq.51.or.itype.eq.45) then
             write (line(ip),751) date0,date1(8:15),isat,ista,cpass(ip),
     |       nobs,iobs,rmean*100,rrms*100,fit*100,a*100,b*1e6,mark
	   else if (itype.eq.40.or.itype.eq.46) then
             write (line(ip),740) date0,date1(8:15),isat,ista,cpass(ip),
     |       nobs,iobs,rmean*1000,rrms*1000,fit*1000,a*1000,b*1e6,mark
	   else if (itype.eq.99.and.altim) then
             write (line(ip),799) date0,date1(8:15),isat,ista,cpass(ip),
     |       nobs,iobs,rmean*100,rrms*100,fit*100,a*100,b*1e6,mark
           else if (itype.eq.100.and.altim) then
             write (line(ip),799) date0,date1(8:15),isat,ista,cpass(ip),
     |       nobs,iobs,rmean*100,rrms*100,fit*100,a*100,b*1e6,mark
	   endif
	 endif
	 l=lnblnk(line(ip))
	 if (.not.long) l=min(l,77)
	 if (l.ne.0) write (*,550) line(ip)(:l)
      enddo

* Process statistics per station/satellite/datatype combination

      write (*,550)
      iobs=0
      nobs=0
      rmean=0
      rrms=0
      kpass=0
      ista=0
      isat=0
      itype=0
      do ii=1,npass
	 ip=isort(ii,1)
	 if (ipass(3,ip).gt.0) then
	    if ((ista.ne.ipass(1,ip).or.isat.ne.ipass(2,ip).or.
     |		 itype.ne.ipass(4,ip)) .and. nobs.ne.0)
     |		 call prstat(iobs,rmean,rrms,itype,kpass,isat,ista,
     |			cpass(iplast),nobs,long,sitename(iplast))
	    rmean=rmean+pobs(ip)*pmean(ip)
	    rrms =rrms +pobs(ip)*prms(ip)**2
	    iobs =iobs +pobs(ip)
	    nobs =nobs +ipass(3,ip)
	    kpass=kpass+1
	    ista=ipass(1,ip)
	    isat=ipass(2,ip)
	    itype=ipass(4,ip)
	    iplast=ip
	 endif
      enddo
      if (nobs.ne.0) call prstat(iobs,rmean,rrms,itype,kpass,isat,ista,
     |			cpass(iplast),nobs,long,sitename(iplast))
      write (*,550)

* Process statistics per satellite / measurement combination

      do i=1,ntypes
	 itype=itypes(i)
	 do j=1,nsats
	    isat=isats(j)
            do ii=1,npass
               ip=isort(ii,1)
               if (ipass(3,ip).gt.0 .and. isat.eq.ipass(2,ip)
     |			.and. itype.eq.ipass(4,ip)) then
	          rmean=rmean+pobs(ip)*pmean(ip)
	          rrms =rrms +pobs(ip)*prms(ip)**2
	          iobs =iobs +pobs(ip)
	          nobs =nobs +ipass(3,ip)
	          kpass=kpass+1
	          if (ista.ne.ipass(1,ip)) nsta=nsta+1
	          ista=ipass(1,ip)
	          iplast=ip
	       endif
	    enddo
            if (nobs.ne.0)
     |call prstat(iobs,rmean,rrms,itype,kpass,isat,nsta,
     |		typenm1(itype),nobs,long,typenm2(itype))
	 enddo
      enddo

*400   continue

* If plot requested....

      if (dev.ne.' ') then

* Make a plot using the PGPLOT routines

      if (pgbeg(0,dev,1,1).ne.1)
     |		call fin("lres: unable to open plot device.")
      call pgask(.false.)

* Set colours

      do i=1,6
	 call pgscr(i-1,rgb(1,i),rgb(2,i),rgb(3,i))
      enddo

      select=.false.
      ii=0
      di=1
410	 ii=ii+di
         di=di/abs(di)
	 if (ii.lt.1) then
	    ii=1
	    di=+1
	    write (*,551) char(7)
	 else if (ii.gt.npass) then
	    ii=npass
	    di=-1
	    write (*,551) char(7)
	 endif
	 ip=ii
	 if (.not.timesort) ip=isort(ii,1)
	 if (line(ip).eq.' ') goto 410
	 l=lnblnk(line(ip))
	 if (select.and.line(ip)(l:l).ne.'!') goto 410
	 ista=ipass(1,ip)
         nobs=ipass(3,ip)
* Mod. by Eelco Doornbos: Only edit range data
*         if(ipass(4,ip).eq.51 .or. ipass(4,ip).eq.45) then
         if(ipass(4,ip).eq.51) then
         call resedit(nobs,resid(1,ip),sigma(1,ip),deriv(1,ip),
     |                elev(1,ip),line(ip),sitename(ip))
*         else if(ipass(4,ip).eq.46.or.ipass(4,ip).eq.40) then
*         call dopresedt(nobs,times(1,ip),resid(1,ip),sigma(1,ip),
*     |       deriv(1,ip),elev(1,ip),line(ip),sitename(ip))
         end if
* End of mod.
	 if (di.eq.0 .or. (curs.eq.0 .and. ii.eq.npass)) goto 500
      goto 410

500   call pgend

      endif	! End of plotting routines

* If delete cards are requested, store them

      if (delfile.ne.' ' .and. returncode.ne.999) then
	 open (40,file=delfile)
	 close (40,status='delete')
	 open (40,file=delfile)
	 do ip=1,npass
	    ista=ipass(1,ip)
            nobs=ipass(3,ip)
	    itype=ipass(4,ip)
	    if (.not.altim .and. (itype.eq.99 .or. itype.eq.100))
     |		goto 510
            if (itype.eq.46 .or. itype.eq.45 .or. itype.eq.40)
     |          goto 510
	    idel0=0
	    idel1=0
	    do i=1,nobs
	       if (sigma(i,ip).gt.0) then
	       else if (idel0.eq.0) then
		  idel0=i
		  idel1=i
	       else if (idel1.eq.i-1) then
		  idel1=i
	       else
		  call prdel(ista,itype,
     |			times(idel0,ip)-t1985,times(idel1,ip)-t1985)
		  idel0=i
		  idel1=i
	       endif
	    enddo
	    if (idel0.ne.0) call prdel(ista,itype,
     |			times(idel0,ip)-t1985,times(idel1,ip)-t1985)
510	    continue
	 enddo
	 close (40)
      endif	! End of delete card generation

* End of program
      
      goto 9999

* Formats

  550 format (a)
  551 format (a,$)
  600 format('usage: lres [options] [ filenames(s) ]'//
     |'where [options] are:'//
     |'  -t       : sort by time'/
     |'  -s       : sort by station (default)'/
     |'  -z       : exclude stations without proper observations'/
     |'  -a       : exclude pass-by-pass altimeter residuals'/
     |'  -l       : print long information per pass'/
     |'  -D       : print statistics of the deleted measurements'/
     |'  -f       : adjust editing using data fit'/
     |'  -d[file] : print delete cards to ''file'' (default: lres.del)'/
     |'  -e[dev]  :',
     |' edit passes using plot device ''dev'' (default: /xs)'/
     |'             (option also sets -t -a -l -d)'//
     |'and'//
     |'filename(s): residual file name(s), default=fort.19')
  740 format (a,' - ',a,2x,i7.7,1x,i4,1x,a8,2i6,2f8.3,2f7.3,f7.1,a)
  751 format (a,' - ',a,2x,i7.7,1x,i4,1x,a8,2i6,2f8.1,2f7.1,f7.1,a)
  799 format (a,' - ',a,2x,i7.7,1x,i4,1x,a8,2i6,2f8.1,2f7.1,f7.1,a)

 9999 end

***********************************************************************
* Subroutine meanrms: Compute mean and RMS of a number of data values

      subroutine meanrms(n,mean,rms)
      integer n
      real*8 mean,rms
      if (n.gt.0) then
         mean=mean/n
         rms=sqrt(rms/n)
      endif
      end

***********************************************************************
* Subroutine prstat: Actually print one statistics line

      subroutine prstat
     |(iobs,rmean,rrms,itype,kpass,isat,ista,stanm,nobs,long,sitename)
      integer iobs,itype,kpass,ista,nobs,isat,l,lnblnk
      logical long
      real*8 rmean,rrms
      character stanm*8,sitename*32

      call meanrms(iobs,rmean,rrms)
      if (itype.eq.40.or.itype.eq.46) write (*,840)
     |	   kpass,isat,ista,stanm,nobs,iobs,rmean*1000,rrms*1000
      if (itype.eq.51.or.itype.eq.45) write (*,851)
     |	   kpass,isat,ista,stanm,nobs,iobs,rmean*100,rrms*100
      if (itype.eq.99) write (*,899)
     |	   kpass,isat,ista,stanm,nobs,iobs,rmean*100,rrms*100
      if (itype.eq.100) write (*,899)
     |	   kpass,isat,ista,stanm,nobs,iobs,rmean*100,rrms*100
      l=lnblnk(sitename)
      if (l.gt.0 .and. long) then
         write (*,550) sitename(:l)
      else
         write (*,'(a)')
      endif
      iobs=0
      nobs=0
      rmean=0
      rrms=0
      kpass=0
      ista=0
  550 format (2x,a)
  840 format (22x,i4,2x,i7.7,1x,i4,1x,a8,2i6,2f8.3,$)
  851 format (22x,i4,2x,i7.7,1x,i4,1x,a8,2i6,2f8.1,$)
  899 format (22x,i4,2x,i7.7,1x,i4,1x,a8,2i6,2f8.1,$)
      end

***********************************************************************
* Subroutine types: allocate array to measurement type

      subroutine types(array,m,n,i)
      integer m,array(m),n,i,j

      do j=1,n
	 if (array(j).eq.i) return
      enddo
      n=n+1
      if (n.gt.m) call fin('lres: too many measurement types or sats')
      do j=n-1,1,-1
	 if (array(j).gt.i) then
	    array(j+1)=array(j)
	 else
	    array(j+1)=i
	    return
	 endif
      enddo
      array(1)=i
      end

***********************************************************************
* Subroutine prdel: print the delete line

      subroutine prdel(i,j,t0,t1)
      character*14 date0,date1
      real*8 t0,t1
      integer i,j

      call strf1985(date0,'%y%m%d%H%M%S',int(t0-0.05))
      call strf1985(date1,'%y%m%d%H%M%S',int(t1+1.05))

      write (40,600) i,j,date0,date1
600   format ('DELETE',i8,i3.3,t41,a,t61,a)
      end
