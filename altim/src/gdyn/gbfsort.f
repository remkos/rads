      program gbfsort

* This program sorts all data in a GBF file timewise, and eliminates double
* occurrances of the same pass (later occurrances overrule earlier).
* For this to work correctly, the data should be either:
* - sorted passwise (not necessarily timewise)
* - timewise (not necessarily passwise)
* - a combination of the above
*
* syntax: gbfsort GBFin ... GBFout
*
* gbfsort can also be called as gbflist, in which case it only gives a pass
* inventory.
*
* syntax: gbflist GBFin ...
*
*-
*  3-Aug-1998 - Created by Remko Scharroo
*  3-Nov-1999 - Added DATEARG function
* 12-Jan-2000 - Added indications 'UNIQ' and 'USED' in output. Fixed
*               proper identification of DBLE passes wit 1 NP.
* 15-Jan-2002 - Changed output of SatID from "i10" to "3x,i7.7"
*-----------------------------------------------------------------------

      include 'gbfsort.inc'
      integer*4 iarg,iargc,narg,ios
      real*8	dum
      character*80 arg
      logical*4 datearg

* Check command line arguments

      call getarg(0,arg)
      if (index(arg,'gbflist').gt.0) then
         mode=0
         narg=iargc()
      else
         mode=1
         narg=iargc()-1
      endif
      if (narg.gt.0) then
      else if (mode.eq.0) then
         write (0,600)
         write (0,602)
	 goto 9999
      else
         write (0,601)
         write (0,602)
	 goto 9999
      endif
600   format('gbflist - List passes in GBF file(s)'//
     |'syntax: gbflist [options] GBFin ...')
601   format('gbfsort - Sort GBF data in file(s) by time and',
     |' eliminates double passes'//
     |'syntax: gbfsort [options] GBFin ... GBFout')
602   format(/'where [options] are:'/
     |' t=t0,t1  : select time interval'/
     |'        ... or use mjd=, doy=, ymd=, sec='/
     |' -u       : only select unique measurements'/
     |' -l4      : list only station numbers'/
     |' -l3      : list pass, station, and global summary (default)'/
     |' -l2      : list station and global summary'/
     |' -l1      : list global summary only'/
     |' -l0      : no list'/
     |' -l       : same as -l0')

* Do the main stuff

      call init
      do iarg=1,narg
	 call getarg(iarg,arg)
	 if (datearg(arg,t0,t1,dum)) then
	 else if (arg(:2).eq.'-u') then
	    uniq=1
	 else if (arg(:2).eq.'-l') then
	    list=0
	    read (arg(3:),*,iostat=ios) list
	 else
	    open (10,file=arg,form='unformatted',status='old')
            call make_inventory
	    close (10)
         endif
      enddo
      call edit_inventory
      if (list.ge.1) call list_inventory
      if (mode.eq.1) then
	 call getarg(narg+1,arg)
	 open (10,file=arg,form='unformatted',status='new')
         call sort_data
	 close (10)
      endif
      if (list.lt.4) write (*,'(/"Normal end of gbfsort")')

* Bye Bye

9999  end

      subroutine init

* Initialise

      include 'gbfsort.inc'
      integer*4 minint4,maxint4,i
      parameter (maxint4=2147483647,minint4=-maxint4-1)

      do i=1,maxpass
         p_t0(i)=maxint4
         p_t1(i)=maxint4
	 p_dble(i)=0
      enddo
      do i=1,9999
         s_t0(i)=maxint4
	 s_t1(i)=minint4
	 s_used(i)=0
	 s_usnr(i)=0
	 s_dble(i)=0
	 s_dbnr(i)=0
      enddo
      t_t0=maxint4
      t_t1=minint4
      t_used=0
      t_usnr=0
      t_dble=0
      t_dbnr=0
      nmeas=0
      npass=0
      t0=minint4
      t1=maxint4
      list=3
      uniq=0
      end

      subroutine make_inventory

* Read the GBF data and make an inventory

      include 'gbfsort.inc'
      character*68 gbf
      integer*4 idsat,istat,mjd,l
      real*8	tfrac,tmeas
      integer*2 mtype
      equivalence (gbf(1:4),idsat),(gbf(5:6),mtype),(gbf(9:12),istat),
     |		  (gbf(17:20),mjd),(gbf(21:28),tfrac)

* Read measurements from UNIT 10

10    read (10,end=100) gbf
      tmeas=((mjd-46066)+tfrac)*86400d0

* Skip measurements outside time interval

      if (tmeas.lt.t0 .or. tmeas.gt.t1) goto 10

* If the measurement does NOT come directly after the end of the last pass
* of the current station, then it is registered as a new pass.

      nmeas=nmeas+1
      if (nmeas.gt.maxmeas) call fin('gbfsort: too many measurements')
      l=s_used(istat)
      if (l.gt.0 .and.
     |	p_mtype(l).eq.mtype .and. p_idsat(l).eq.idsat .and.
     |	tmeas.gt.p_t1(l) .and. tmeas.le.p_t1(l)+dt) then
	 p_t0(l)=min(p_t0(l),tmeas)
	 p_t1(l)=max(p_t1(l),tmeas)
	 p_nr(l)=p_nr(l)+1
      else
	 npass=npass+1
         if (npass.gt.maxpass) call fin('gbfsort: too many passes')
	 l=npass
	 s_used(istat)=l
	 p_istat(l)=istat
	 p_mtype(l)=mtype
	 p_idsat(l)=idsat
	 p_t0(l)=tmeas
	 p_t1(l)=tmeas
	 p_nr(l)=1
      endif
      m_pass(nmeas)=l
      m_time(nmeas)=tmeas
      m_gbf (nmeas)=gbf
      goto 10
100   return
      end

      subroutine edit_inventory

* Flag the double passes
* A double pass is one with overlapping time span. The earlier one in
* the data set is marked 'DBLE', the later is kept and marked 'USED'.

      include 'gbfsort.inc'
      integer*4 i,j

      do i=npass,2,-1
         if (p_dble(i).ne.1) then
	    do j=1,i-1
	       if (p_istat(j).eq.p_istat(i) .and.
     |		p_idsat(j).eq.p_idsat(i) .and. p_mtype(j).eq.p_mtype(i)
     |		.and. p_t0(j).le.p_t1(i)+dt .and. p_t1(j).ge.p_t0(i)-dt) then
		  p_dble(j)=1; p_dble(i)=-1
	       endif
	    enddo
	 endif
      enddo
      end

      subroutine list_inventory

* Sort the passes by time and then list all passes

      include 'gbfsort.inc'
      character date0*15,date1*15,dble*4
      integer*4 i,j,l

* Initialise

      do i=1,npass
         pnt(i)=i
      enddo
      do i=1,9999
         s_used(i)=0
      enddo

* Sort passes by start time

      call dqsort(pnt,p_t0,npass)

* Print out pass statistics and keep station statistics

      if (list.eq.3) write (*,605)

      do i=1,npass
         j=pnt(i)
	 l=p_istat(j)

	 s_t0(l)=min(s_t0(l),p_t0(j))
	 s_t1(l)=max(s_t1(l),p_t1(j))
	 if (p_dble(j).eq.1) then
	    dble='DBLE'
	    s_dble(l)=s_dble(l)+1
            s_dbnr(l)=s_dbnr(l)+p_nr(j)
	 else if (p_dble(j).eq.0) then
	    dble='UNIQ'
	    s_used(l)=s_used(l)+1
            s_usnr(l)=s_usnr(l)+p_nr(j)
	 else if (uniq.eq.0) then
	    dble='USED'
	    s_used(l)=s_used(l)+1
            s_usnr(l)=s_usnr(l)+p_nr(j)
	 else
	    dble='NOTU'
	    s_dble(l)=s_dble(l)+1
            s_dbnr(l)=s_dbnr(l)+p_nr(j)
	 endif

	 call strf1985(date0,'%y%m%d %H:%M:%S',nint(p_t0(j)))
	 call strf1985(date1,'%y%m%d %H:%M:%S',nint(p_t1(j)))
	 if (list.eq.3)
     |		write (*,610) j,date0,date1,l,p_idsat(j),p_mtype(j),
     |		p_nr(j),dble
      enddo

* Write only station numbers

      if (list.eq.4) then
         do i=1,9999
	    if (s_used(i).gt.0) write (*,'(i4.4)') i
	 enddo
         return
      endif

* Write station statistics and keep global statistics
      
      if (list.ge.2) write (*,615)

      j=0
      do i=1,9999
         if (s_usnr(i).gt.0) then
	    j=j+1
	    t_used=t_used+s_used(i)
	    t_usnr=t_usnr+s_usnr(i)
	    t_dble=t_dble+s_dble(i)
	    t_dbnr=t_dbnr+s_dbnr(i)
	    t_t0  =min(t_t0,s_t0(i))
	    t_t1  =max(t_t1,s_t1(i))

	    call strf1985(date0,'%y%m%d %H:%M:%S',nint(s_t0(i)))
	    call strf1985(date1,'%y%m%d %H:%M:%S',nint(s_t1(i)))
	    if (list.ge.2)
     |		write (*,620) j,date0,date1,i,s_used(i),s_usnr(i),
     |		s_dble(i),s_dbnr(i)
         endif
      enddo

* Print global statistics
	    
      if (list.ge.1) then
         write (*,625)
	 call strf1985(date0,'%y%m%d %H:%M:%S',nint(t_t0))
	 call strf1985(date1,'%y%m%d %H:%M:%S',nint(t_t1))
         if (j.gt.0)
     |		write (*,630) j,date0,date1,t_used,t_usnr,t_dble,t_dbnr
      endif

605   format (/'*** Pass Statistics ***'/
     |'Passnr            Start               End  Stat',
     |'     SatID Type #NP')
610   format (i6,2x,a15,' - ',a15,i6,3x,i7.7,i5,i4,2x,a4)
615   format (/'*** Station Statistics ***'/
     |'     #            First              Last  Stat',
     |'  #PS(USED)#NP  #PS(DBLE)#NP')
620   format (i6,2x,a15,' - ',a15,i6,2(i6,i8))
625   format (/'*** Global Statistics ***'/
     |' #Stat            First              Last      ',
     |'  #PS(USED)#NP  #PS(DBLE)#NP')
630   format (i6,2x,a15,' - ',a15,6x,2(i6,i8))
      end

      subroutine sort_data

* Sort the remaining GBF data by time and write them to file

      include 'gbfsort.inc'
      integer*4 i,n,j

* Initialize

      n=0
      do i=1,nmeas
	 j=p_dble(m_pass(i))
         if (j.ne.1 .and. uniq*j.eq.0) then
	    n=n+1
	    pnt(n)=i
 	 endif
      enddo
      if (n.ne.t_usnr) write (*,*) 'sort: n,t_usnr=',n,t_usnr

* Sort the remaining data

      call dqsort(pnt,m_time,n)
         
* Write GBF data

      do i=1,n
	 write (10) m_gbf(pnt(i))
      enddo
      end
