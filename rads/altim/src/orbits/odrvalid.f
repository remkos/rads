**ODRVALID -- Program to validate ODRs
*
*-
* 30-Oct-1998 - Created
*-----------------------------------------------------------------------
      program	lodr
      integer	i,iarg,iargc,fd,ios,day,mday,mman,nman
      parameter (mday=40,mman=100)
      integer	rep,begin,end,time0,time1,tstep,ver,arc,nrec
      integer	rec(4),nr(mday),yymmdd,mdate,man(mman)
      character odrnm*80,orbdir*80,dev*80,satel*8
      real*8	alt1,alt2,lat,lon,rms(mday),mean(mday),rev
      real	x(mday),y(mday),daymax
      logical	first/.true./

* Scan arguments

      call getarg(1,orbdir)
      call getarg(2,dev)

* Read manoeuvre file

      open (10,file='/user/remko/ers/ODR.ERS-2/maneuver.txt')
      do i=1,14
         read (10,*)
      enddo
      nman=0
10    read (10,'(5x,i6)',end=20) yymmdd
      nman=nman+1
      man(nman)=(mdate(2,yymmdd)-46066)*86400
      write(*,*) nman,man(nman)
      goto 10
20    continue

* Open ODR and get header information

      do iarg=3,iargc()
      call getarg(iarg,odrnm)
      call odrinfo(fd,odrnm,satel,rep,arc,ver,nrec,
     |time0,time1,tstep,begin,end,rev)

* Check manoeuvre

      do i=1,nman
         if (man(i).ge.time0 .and. man(i).le.time1) then
	    write (*,*) 'skipped ',odrnm
	    goto 100
         endif
      enddo

* Initialise statistics

      do i=1,mday
	 mean(i)=0
         rms(i)=0
	 nr(i)=0
      enddo

* Process data records

      do i=1,nrec
         day=(i-1)/1440+1
         read (fd,rec=i+2,iostat=ios) rec
         if (ios.ne.0) call fin("premature end of file")
	 alt1=rec(4)*1d-3
         call getorb(rec(1)*1d0,lat,lon,alt2,orbdir,.true.)
	 rms(day)=rms(day)+(alt1-alt2)**2
	 mean(day)=mean(day)+(alt1-alt2)
	 nr(day)=nr(day)+1
      enddo
      do day=1,20
         if (nr(day).gt.1) then
	    write (*,*)	day,nr(day),sqrt(rms(day)/nr(day))*1d2
	    x(day)=day-0.5
	    y(day)=sqrt(rms(day)/nr(day))
            daymax=day
         endif
      enddo

* Start plot

      if (first) then
         first=.false.
         call pgbeg(0,dev,1,1)
	 call pgvport(0.10,0.95,0.10,0.95)
         call pgswin(0.0,daymax,0.0,100.0)
	 call pgsci(15)
	 call pgrect(6.0,7.0,0.0,100.0)
	 call pgsci(1)
         call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
	 call pglab('time from start of orbital arc (days)',
     |		'RMS difference between operational and precise orbit (cm)',
     |		' ')
         call pgswin(0.0,daymax,0.0,1.0)
	 call pgslw(3)
	 call pgmtext('t',1.0,0.1,0.0,'data period')
	 call pgmtext('t',1.0,0.5,0.5,'no SLR data')
	 call pgmtext('t',1.0,0.9,1.0,'prediction')
	 call pgslw(2)
      endif
      call pgline(nint(daymax),x,y)
100   continue
      enddo
      call pgend
      goto 9999

9999  end
