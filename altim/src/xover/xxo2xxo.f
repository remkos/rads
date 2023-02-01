      program xxo2xxo

* Manipulate XXO files. Select and/or add new orbits
*-
* $Log: xxo2xxo.f,v $
* Revision 1.6  2011/02/10 18:40:55  rads
* - Allow up to 5000000 xovers
*
* (c) Remko Scharroo - TU Delft, Altimetrics
*----------------------------------------------------------------------
      implicit  none
      real*8	rmslon/0d0/,rmslat/0d0/,rmsh/0d0/,rmshpr/0d0/
      real*8	rlon0,rlat0
      real*8	rlat1,rlon1,dssh1,dorb1,orbit1,utc1
      real*8	rlat2,rlon2,dssh2,dorb2,orbit2,utc2
      real*8	dellon/0d0/,dellat/0d0/,delh/0d0/,t0,t1,dum
      real*8	pi,rad,latfact,coslat,tb
      real*8	dhellips
      character file1*80/' '/,file2*80/' '/,arg*80,spec*4,
     |		dir*80/' '/,dev*80/' '/
      integer*2 i2(10),it(2)/2*0/,dlt(2)/2*0/
      integer*4 i4(6),ssh(2),orb(2),getorb,tbias/0d0/,npar,maxtrk,
     |		sat1/0/,sat2/0/,trk0/0/,trk1/0/,trk2/0/,trk3/0/,
     |          lon0,lon1,lat0,lat1
      parameter (maxtrk=32767)
      integer*2 seltrk(maxtrk)
      integer*4 maxx,l,it0,it1,lnblnk
      parameter (maxx=5000000)
      real*4	x(maxx,2),o(maxx,2),t(maxx,2),y(maxx,2)
      real*8	tag(maxx*2),rlat(maxx*2),rlon(maxx*2)
      integer*4	key(maxx*2),pnt(maxx*2)
      real*4	xmax/-1e35/,xmin/1e35/,tmax/-1e35/,tmin/1e35/
      logical	ers/.true./,tp/.false./,noreplace/.false./,datearg,
     |		verbose/.false./,high/.false./,low/.false./,ellips/.false./,
     |		change,usedelta/.false./,useaux/.false./
      integer*4	iargc,narg,iarg,irec,nrec,mrec,i,j,
     |		ierr/0/,npnt,mode/0/
      equivalence (arg,it)
      integer*4  minint4,maxint4
      parameter  (maxint4=2147483647,minint4=-maxint4-1)

* Initialize

      lon0=minint4
      lon1=maxint4
      lat0=minint4
      lat1=maxint4
      it0=minint4
      it1=maxint4

* Scan arguments

      narg=iargc()
      if (narg.lt.3) goto 1300
      do iarg=1,narg
	 call getarg(iarg,arg)
	 if (arg(:4).eq.'dev=') then
	    dev=arg(5:)
	 else if (arg(:2).eq.'-a') then
	    ierr=-2
	 else if (arg(:2).eq.'-A') then
	    ierr=-1
	 else if (arg(:2).eq.'-n') then
	    noreplace=.true.
	 else if (arg(:2).eq.'-v') then
	    verbose=.true.
	 else if (arg(:2).eq.'-h') then
	    high=.true.
	 else if (arg(:2).eq.'-l') then
	    low=.true.
	 else if (arg(:2).eq.'-d') then
	    mode=2
	 else if (arg(:4).eq.'orb=') then
	    dir=arg(5:)
	 else if (arg(:4).eq.'Orb=') then
	    dir=arg(5:)
	    ellips=.true.
	 else if (datearg(arg,t0,t1,dum)) then
	    it0=nint(t0)
	    it1=nint(t1)
	 else if (arg(:4).eq.'lon=') then
	    read (arg(5:),*) rlon0,rlon1
	    lon0=nint(rlon0*1d6)
	    lon1=nint(rlon1*1d6)
	 else if (arg(:4).eq.'lat=') then
	    read (arg(5:),*) rlat0,rlat1
	    lat0=nint(rlat0*1d6)
	    lat1=nint(rlat1*1d6)
	 else if (arg(:4).eq.'trk=') then
	    read (arg(5:),*,iostat=i) trk0,trk1,trk2,trk3
	    if (trk2.eq.0) trk2=trk0
	    if (trk3.eq.0) trk3=trk1
	 else if (arg(:4).eq.'sat=') then
	    read (arg(5:),*,iostat=i) sat1,sat2
	    if (sat2.eq.0) sat2=sat1
	 else if (arg(:6).eq.'tbias=') then
	    read (arg(7:),*) tb
	    tbias=nint(tb*1d3)
	    write (*,*) "tbias=",tbias
	 else if (file1.eq.' ') then
	    file1=arg
	 else if (file2.eq.' ') then
	    file2=arg
	 else
	    dir=arg
	 endif
      enddo

* Check if sufficient file names or options are supplied.

      if (file1.eq.' ' .or. file2.eq.' ') goto 1300
      if (dir.eq.' ') then
         if (it0.lt.0 .and. trk0.le.0 .and. sat1.le.0) goto 1300
	 mode=1
      endif

* Set some parameters

      pi=4*atan(1d0)
      rad=pi/180
      latfact=rad*6371d3

* If selection on satellite number is requested, XTF file has to be loaded.

      if (sat1.ne.0) then
         l=index(file1,'.xxo')
	 arg=file1
	 arg(l:)='.xtf'
	 open (10,file=arg,status='old',form='unformatted',
     |		access='direct',recl=12)
         read (10,rec=1,err=1305) spec,nrec,npar
	 if (spec(:3).ne.'@XT') call fin('XTF file not in proper format')
	 if (nrec.gt.maxtrk) call fin('Too many tracks')
	 close (10)
	 l=34+npar*8
	 open (10,file=arg,status='old',form='unformatted',
     |		access='direct',recl=l)
         do i=2,nrec+1
	    read (10,rec=i,err=1305) it
	    if (it(2).eq.sat1) seltrk(it(1))=1
	    if (it(2).eq.sat2) seltrk(it(1))=3-seltrk(it(1))
	 enddo
	 close (10)

* If no satellite is selected, but tracknumbers, select those.

      else if (trk0.ne.0) then
	 do j=trk0,trk1
	    seltrk(j)=1
	 enddo
	 do j=trk2,trk3
	    seltrk(j)=3-seltrk(j)
	 enddo

* If nothing is selected, select all.

      else
	 do j=1,maxtrk
            seltrk(j)=2
	 enddo
      endif

* Open input XXO file

12    open (10,file=file1,recl=44,access='direct',status='old',
     |		form='unformatted')
      read (10,rec=1) spec,nrec
      if (spec.eq.'@D_H') then
	 close (10)
         open (11,file=file1,recl=4,access='direct',status='old',
     |		form='unformatted')
	 l=index(file1,'.xxo')
	 file1=file1(:l+3)
	 usedelta=.true.
	 goto 12
      else if (spec.ne.'@XXO') then
	 call fin('xxo2xxo: input file not XXO')
      endif

* Open output delta file (mode=2) or XXO file (mode<>2)

      if (mode.eq.2) then
         open (20,file=file2,recl=4,access='direct',
     |	status='new',form='unformatted')
      else
         open (20,file=file2,recl=44,access='direct',
     |	status='new',form='unformatted')
      endif

* If merely selection on time, tracknr, or satellite is requested:

      if (mode.eq.1) then

* Check if AUX file exists, open input and output AUX file.

      l=lnblnk(file1)
      arg=file1
      arg(l+1:)='.aux'
      inquire (exist=useaux,file=arg)
      if (useaux) then
         open (12,file=arg,recl=12,status='old',form='unformatted',
     |		access='direct')
         l=lnblnk(file2)
	 npar=6
         arg=file2
         arg(l+1:)='.aux'
         open (22,file=arg,recl=12,status='new',form='unformatted',
     |		access='direct')
      endif
     
      mrec=0
      do irec=1,nrec
	 read (10,rec=irec+1,err=1305) i4,it,ssh,orb

* Here we use a trick:
* - If one satellite is selected, all selected track numbers i have
* seltrk(i)=2, so when the sum is 4. Combinations of selected and non
* selected tracks are 2 or 0.
* - If two satellites (or two batches of tracks) are selected, one satellite
* (or one batch) gets seltrk(i)=1, and the other seltrk(i)=3. Now only the
* cross-combination of the two satellites or batches get a sum of 4.
* - If satellites are not specified, all seltrk(i)=2, so the sum is always 4.

	 if (lat0.le.i4(1) .and. i4(1).le.lat1 .and.
     |	     lon0.le.i4(2) .and. i4(2).le.lon1 .and.
     |	     it0.le.i4(3) .and. i4(3).le.it1 .and.
     |	     it0.le.i4(5) .and. i4(5).le.it1 .and.
     |	    seltrk(it(1))+seltrk(it(2)).eq.4) then

* Add delta to orbital altitude and ssh if requested

	    if (usedelta) then
	       read (11,rec=irec+6,err=1305) dlt
	       do j=1,2
	          ssh(j)=ssh(j)+dlt(j)*1000
	          orb(j)=orb(j)+dlt(j)
	       enddo
	    endif
	    mrec=mrec+1
	    write (20,rec=mrec+1) i4,it,ssh,orb

* If AUX file available, select also this one

	    if (useaux) then
	       read (12,rec=irec+1) (i2(j),j=1,npar)
	       write (22,rec=mrec+1) (i2(j),j=1,npar)
	    endif
	 endif
      enddo
      write (20,rec=1) spec,mrec
      write (*,610) mrec,nrec
      if (useaux) then
         write (22,rec=1) '@AUX',mrec,npar
      endif
610   format ('xxo2xxo: ',i9,' records selected of ',i9)
      goto 9999
      endif	! (mode.ne.0)
	 
      if (index(dir,'TOPEX').gt.0 .or. index(dir,'JASON').gt.0) then
	 tp=.true.
	 ers=.false.
      endif

* First store the time tags and sort them
* Apply time tag bias

      npnt=0
      do irec=1,nrec
	 read (10,rec=irec+1,err=1305) i4,it,ssh,orb
	 i4(4)=i4(4)-tbias
	 i4(6)=i4(6)-tbias
	 if (change(high,low,it(1).gt.it(2),orb(1),ers,tp)
     |	.and. i4(3).ge.it0 .and. i4(3).le.it1) then
	    npnt=npnt+1
            if (npnt.gt.maxx*2) goto 1320
	    pnt(npnt)=npnt
	    key(npnt)=i4(3)
	    tag(npnt)=i4(3)+i4(4)/1d6
	 endif
	 if (change(high,low,it(2).gt.it(1),orb(2),ers,tp)
     |	.and. i4(5).ge.it0 .and. i4(5).le.it1) then
	    npnt=npnt+1
            if (npnt.gt.maxx*2) goto 1320
	    pnt(npnt)=npnt
	    key(npnt)=i4(5)
	    tag(npnt)=i4(5)+i4(6)/1d6
	 endif
      enddo
      if (verbose) write (0,551) 'Sorting time tags ...'
      call iqsort(pnt,key,npnt)
      if (verbose) write (0,550) ' Done'

* Compute new orbit for each time tag

      do irec=1,npnt
	 j=pnt(irec)
	 i=getorb(tag(j),rlat(j),rlon(j),tag(j),dir,verbose)
	 if (i.gt.0 .or. i.lt.ierr) then
	    tag(j)=0
	 else if (ellips) then
	    tag(j)=tag(j)+dhellips(3,rlat(j))
         endif
      enddo

* Reread all data files and replace orbit

      npnt=0
      mrec=0
      do irec=1,nrec

* Read XXO record. Read D_H record if required.

*	 write (*,*) irec,nrec
	 read (10,rec=irec+1,err=1305) i4,it,ssh,orb
	 if (usedelta) then
	    read (11,rec=irec+6,err=1305) dlt
	 else
	    dlt(1)=0
	    dlt(2)=0
	 endif
	 utc1=i4(3)+i4(4)/1d6
	 utc2=i4(5)+i4(6)/1d6

	 if (change(high,low,it(1).gt.it(2),orb(1),ers,tp)
     |	.and. i4(3).ge.it0 .and. i4(3).le.it1) then
	    npnt=npnt+1
	    if (noreplace .and. tag(npnt).ne.0) then
	       orbit1=(orb(1)+dlt(1))/1d3
	       rlat1=i4(1)/1d6
	       rlon1=i4(2)/1d6
	    else
	       orbit1=tag(npnt)
	       rlat1=rlat(npnt)
	       rlon1=rlon(npnt)
	       i4(4)=i4(4)-tbias
	    endif
	 else
	    orbit1=(orb(1)+dlt(1))/1d3
	    rlat1=i4(1)/1d6
	    rlon1=i4(2)/1d6
	 endif

	 if (change(high,low,it(2).gt.it(1),orb(2),ers,tp)
     |	.and. i4(5).ge.it0 .and. i4(5).le.it1) then
	    npnt=npnt+1
	    if (noreplace .and. tag(npnt).ne.0) then
	       orbit2=(orb(2)+dlt(2))/1d3
	       rlat2=i4(1)/1d6
	       rlon2=i4(2)/1d6
	    else
	       orbit2=tag(npnt)
	       rlat2=rlat(npnt)
	       rlon2=rlon(npnt)
	       i4(6)=i4(6)-tbias
	    endif
	 else
	    orbit2=(orb(2)+dlt(2))/1d3
	    rlat2=i4(1)/1d6
	    rlon2=i4(2)/1d6
	 endif

	 if (orbit1.eq.0 .or. orbit2.eq.0) then
	    if (mode.eq.2) then
	       dlt(1)=-32768
	       dlt(2)=-32768
	       write (20,rec=irec+6) dlt
	    endif
	 else

	    dorb1=orbit1-orb(1)/1d3
	    dssh1=ssh(1)/1d6+dorb1
	    dorb2=orbit2-orb(2)/1d3
	    dssh2=ssh(2)/1d6+dorb2

	    mrec=mrec+1
	    t(mrec,1)=utc1/86400.
	    t(mrec,2)=utc2/86400.
	    x(mrec,1)=(ssh(1)-ssh(2))/1d4
	    x(mrec,2)=-x(mrec,1)
	    y(mrec,1)=(dssh1-dssh2)*100
	    y(mrec,2)=-y(mrec,1)
	    o(mrec,1)=dorb1*100
	    o(mrec,2)=dorb2*100
	    tmin=min(t(mrec,1),t(mrec,2),tmin)
	    tmax=max(t(mrec,1),t(mrec,2),tmax)
	    xmin=min(xmin,x(mrec,1),x(mrec,2),o(mrec,1),o(mrec,2),
     |		y(mrec,1),y(mrec,2))
	    xmax=max(xmax,x(mrec,1),x(mrec,2),o(mrec,1),o(mrec,2),
     |		y(mrec,1),y(mrec,2))
	    coslat=cos(rlat2*rad)
*	 write (*,'(6f12.6)') rlat1,rlat2,(rlat1-rlat2)*latfact,
*     .rlon1,rlon2,(rlon1-rlon2)*latfact*coslat
	    rmslon=rmslon+((rlon1-rlon2)*coslat)**2
	    rmslat=rmslat+((rlat1-rlat2)       )**2
	    rmsh  =rmsh  +(dssh1-dssh2)**2
	    rmshpr=rmshpr+(ssh(1)/1d6-ssh(2)/1d6)**2
            rlon2=(rlon1+rlon2)/2
            rlat2=(rlat1+rlat2)/2
	    dellon=dellon+((i4(2)/1d6-rlon2)*coslat)**2
	    dellat=dellat+((i4(1)/1d6-rlat2)       )**2
	    delh  =delh  +dorb1**2+dorb2**2
	    if (mode.eq.0) then
	       orb(1)=nint(orbit1*1d3)
	       orb(2)=nint(orbit2*1d3)
	       ssh(1)=nint(dssh1*1d6)
	       ssh(2)=nint(dssh2*1d6)
	       i4(1)=nint(rlat2*1d6)
	       i4(2)=nint(rlon2*1d6)
	       write (20,rec=mrec+1,err=1310) i4,it,ssh,orb
	    else
	       dlt(1)=nint(dorb1*1d3)
	       dlt(2)=nint(dorb2*1d3)
	       write (20,rec=irec+6,err=1310) dlt
	    endif
         endif
      enddo
      close (10)
      if (mode.eq.2) then
	 write (20,rec=1) '@D_H'
	 write (20,rec=2) spec
	 write (20,rec=3) nrec
	 it(1)=2
	 it(2)=0
	 write (20,rec=4) it
	 it(1)=9
	 it(2)=11
	 write (20,rec=5) it
	 it(1)=10
	 it(2)=12
	 write (20,rec=6) it
      else
         write (20,rec=1) spec,mrec
      endif
      close (20)
      write (*,600) nrec,mrec,
     |sqrt(dellon/mrec)*latfact,sqrt(dellat/mrec)*latfact,
     |sqrt(delh/mrec/2),
     |sqrt(rmslon/mrec)*latfact,sqrt(rmslat/mrec)*latfact,
     |sqrt(rmsh/mrec),sqrt(rmshpr/mrec)

600   format (
     |'recs in = ',i7/
     |'recs out= ',i7/
     |'---------- CHANGES TO XOVER FILE -----'/
     |'lon RMS = ',f10.6,' m'/
     |'lat RMS = ',f10.6,' m'/
     |'h   RMS = ',f10.6,' m'/
     |'---------- XOVER DIFFERENCES ---------'/
     |'lon RMS = ',f10.6,' m'/
     |'lat RMS = ',f10.6,' m'/
     |'h   RMS = ',f10.6,' (',f10.6,') m')

      if (dev.eq.' ') goto 9999
      call pgbeg(0,dev,1,1)
      call pgswin(tmin,tmax,xmin,xmax)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      i=index(file1,' ')-1
      call pglabel('days (1985)','cm',file1(:i)//' with '//
     .	dir)
      call pgpt(mrec,t(1,1),o(1,1),2)
      call pgpt(mrec,t(1,2),o(1,2),2)
      call pgsci(2)
      call pgpt(mrec,t(1,1),x(1,1),5)
      call pgpt(mrec,t(1,2),x(1,2),5)
      call pgsci(3)
      call pgpt(mrec,t(1,1),y(1,1),5)
      call pgpt(mrec,t(1,2),y(1,2),5)
      call pgend
      goto 9999

  550 format (a)
  551 format (a,$)
 1300 write (0,1301)
 1301 format ('xxo2xxo rewrites xxo-files'//
     |'usage: xxo2xxo [options] XXO-file-in XXO-file-out'//
     |'where [options] are'/
     |'dev=dev    : make a plot on device "dev"'/
     |'orb=orbdir : specify orbit directory'/
     |'Orb=orbdir : specify orbit directory (+ ellips conv.)'/
     |'lon=x0,x1  : specify longitude boundaries (deg)'/
     |'lat=y0,y1  : specify latitude  boundaries (deg)'/
     |'t=t0,t1    : select time period ([yy]yymmdd[hhmmss],mjd)'/
     |'         ... or use mjd=, doy=, ymd=, sec='/
     |'sat=sat    : select satellite sat'/
     |'sat=sat1,sat2 : select dual-satellite xovers'/
     |'-v         : give some processing info'/
     |'-arc       : use the entire arc'/
     |'-Arc       : use only first 5.5 days of arc'/
     |'-noreplace : do not replace orbit'/
     |'-low       : replace only orbit on lower track number'/
     |'-high      : replace only orbit on higher track number'/
     |'-d         : create a @D_H (delta) file iso full XXO'/
     |'tbias=dt   : subtract time tag bias dt (msec)')
      goto 9999
 1305 call fin('xxo2xxo: error reading input file')
 1310 call fin('xxo2xxo: error writing output file')
 1320 call fin('xxo2xxo: too many input records')
 9999 end

      function change(high,low,higher,orb,ers,tp)
      logical change,high,low,higher,ers,tp
      integer*4 orb,thkm
      parameter (thkm=1 000 000 000)
      if (high) then
	 change=higher
      else if (low) then
	 change=.not.higher
      else if (orb.lt.thkm.and.ers) then
	 change=.true.
      else if (orb.gt.thkm.and.tp) then
	 change=.true.
      else
	 change=.false.
      endif
      end
