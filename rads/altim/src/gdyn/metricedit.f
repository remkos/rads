**METRICEDIT -- Edit a METRIC file
*
* METRICEDIT is a program to select a part of a METRIC file 
* based on time limits. Also, it can set the station
* depended SIGMAS (weights) to the appropriate values and excluded
* measurements from the orbit determination process by setting the
* SIGMA to a negative value.
*-
* 27-Mar-2002 -- Created by Eelco Doornbos from GBFEDIT 
*-----------------------------------------------------------------------
      program metricedit
      implicit none

* Arguments

      integer maxfiles,nfile/0/,ifile
      parameter (maxfiles=50)
      character*80 arg,delfile/' '/,sigfile/' '/,filenm(maxfiles)
      character*15 date1,date2,date3,date4
      integer iargc,iarg,ios,combi(10)/10*0/,mcombi
      real*8 t0/-1d30/,t1/1d30/,setrange/-1d0/
      logical datearg

* Sigma and delete cards

      integer del0,del1,nread,l
      real*4 sigma0,sigma(10000,100)
      include 'gbfedit.inc'

* Station statistics

      logical stat/.false./,pass/.false./,warn(10000)
      integer s_pnt(10000),s_del(10000),s_pas(10000),s_pts(10000)
      real*8  s_last(10000),s_first(10000)
      integer ncombi(10),lcombi/0/
      real*8  tcombi/0d0/

* METRIC variables

      integer*4 isat1, istat1, isat2, istat2, isat3, istat3
      integer*4 ityptim, imet1,imet2,imet3, ipre
      integer*4 itime
      real*8    tobs1, tobs2, range, interval
      real*8    tdrs1,tdrs2,tdrs3,tdrs4
      real*8    antpre
      real*8    cmass,axis,ctrop,cion,ctrans,crel,dum1
      real*8    vlight,wavel,sig
  
      integer*4 mtype

      real*8 t,tmin1,tmax1,tmin2,tmax2

* General

      integer i,j,nsel
      real*8 tpass,dum(4)
      parameter (tpass=1000/86400d0)

* Scan arguments

      nfile=0
      do iarg=1,iargc()
         call getarg(iarg,arg)
         if (arg(:4).eq.'del=') then
            delfile=arg(5:)
         else if (arg(:4).eq.'sig=') then
            sigfile=arg(5:)
	 else if (arg(:5).eq.'sig0=') then
	    read (arg(6:),*) sigma0   
	    sigma0=sigma0/100
	 else if (datearg(arg,t0,t1,dum)) then
            t0=t0/86400+46066
            t1=t1/86400+46066
	 else if (arg(:6).eq.'range=') then
	    read (arg(7:),*) setrange
	 else if (arg(:2).eq.'-p') then
	    pass=.true.
	 else if (arg(:2).eq.'-s') then
	    stat=.true.
	 else if (arg(:6).eq.'combi=') then
	    read (arg(7:),*,iostat=ios) combi
	    do i=1,10
	       if (combi(i).gt.0) mcombi=i
	    enddo
         else
            nfile=nfile+1
	    if (nfile.gt.maxfiles)
     |		call fin("metricedit: too many input files!")
            filenm(nfile)=arg
         endif
      enddo
      if (nfile.lt.2) then
         write (0,10)
         goto 9999
      endif
10    format ('METRICEDIT -- Edit METRIC files'//
     |'syntax: metricedit [options] infile(s) outfile'//
     |'where'//
     |'infile(s) : input METRIC file(s)'/
     |'outfile   : output METRIC file'//
     |'and [options] are:'/
     |'t=t0,t1   : select period t0,t1 (yymmdd,yyddd)'/
     |'        ... or use mjd=, doy=, ymd=, sec='/
     |'del=file  : read DELETE cards from "file"'/
     |'sig=file  : read SIGMA cards from "file"'/
     |12x,'(defaults will be read from the system.data file)'/
     |'sig0=sig0 : add (rss) "sig0" to all SIGMA (cm)'/
     |'range=rng : set ranges to "rng" (m)'/
     |'-stat     : write statistics per station')

* Initialise

      do i=1,10000
       do j=1,100
         sigma(i,j)=1
	 pnt1(i)=0
	 warn(i)=.false.
       enddo
      enddo
       
* Read sigma cards

      if (sigfile.ne.' ') then
         open (10,file=sigfile,status='old')
200      read (10,'(a)',end=299) arg
         read (arg,220) i,j,sigma(i,j)
220      format (10x,i4,1x,i2,13x,f10.8)
         goto 200
299      close (10)
      endif

* Read delete file

      ndel=0
      if (delfile.ne.' ') call readdel(delfile)    

* Open output METRIC

      open (20,file=filenm(nfile),form='unformatted')

* Read and edit METRIC files

      nsel=0
      ndel=0
      del0=1
      del1=0
      tmin2=1d40
      tmax2=-1d40
      do ifile=1,nfile-1
         open (10,file=filenm(ifile),status='old',form='unformatted')
         nread=0
         tmin1=1d40
         tmax1=-1d40
400      read (10,end=499)
     |    isat1, istat1, isat2, istat2, isat3, istat3,
     |    ityptim, imet1,imet2,imet3, ipre,
     |    itime, tobs1, tobs2, range, interval,
     |    tdrs1,tdrs2,tdrs3,tdrs4, antpre,
     |    cmass,axis,ctrop,cion,ctrans,crel,dum1,
     |    vlight,wavel,sig
         print *, sig
         mtype = int(ityptim/1d6)
         t = itime/86400.0D0 + 39855.0D0 + tobs1
	 tmin1=min(tmin1,t)
	 tmax1=max(tmax1,t)
	 nread=nread+1
         if (t.lt.t0 .or. t.gt.t1) goto 400
	 tmin2=min(tmin2,t)
	 tmax2=max(tmax2,t)
         nsel=nsel+1

* Assign sigma: This comes either from the sigma file (sigma(istat1)<100cm)
* or from the system.data file.
*
* If the original sigma is negative, set the output sigma also to negative.
*
* If the inserted |sigma| > 99cm, give a warning

*         if (sigma(istat1,mtype).lt..9999) then
*	    s=sqrt(sigma(istat1,mtype)**2+sigma0**2)
*	 else
*	    call statinfo(t,istat1,dum,dum,dum,i,dum,dum,dum)
*	    s=sqrt((i/1e2)**2+sigma0**2)
*	 endif 
*	 if (abs(s).gt..980) warn(istat1)=.true.
*	 if (sig.lt.0 .or. s.lt.0) then
*	    sig=-abs(s)
*	 else
*	    sig=s
*	 endif

* Check if the data was marked 'DELETE'

         do i=pnt2(istat1),pnt1(istat1+1)-1
	    pnt2(istat1)=i
	    if (t.lt.delt0(i)) then
	       goto 430			! We're not yet there
	    else if (t.le.delt1(i)) then
	       sig=-abs(sig)		! Remove measurement
	       goto 430
	    else
		! Proceed one DELETE card   
	    endif   
	 enddo
430      continue
	 
* Set range to requested value

	 if (setrange.ge.0) range=setrange
	 
* Write METRIC line to output file

         write(20) 
     |    isat1, istat1, isat2, istat2, isat3, istat3,
     |    ityptim, imet1,imet2,imet3, ipre,
     |    itime, tobs1, tobs2, range, interval,
     |    tdrs1,tdrs2,tdrs3,tdrs4, antpre,
     |    cmass,axis,ctrop,cion,ctrans,crel,dum1,
     |    vlight,wavel,sig

* Update statistics

         if (sig.lt.0) then
	    ndel=ndel+1
	    s_del(istat1)=s_del(istat1)+1
	 else
	    s_pnt(istat1)=s_pnt(istat1)+1
	    if (t-s_last(istat1).gt.tpass) then
	       if (pass)
     |		call writepass(istat1,s_first(istat1),
     |                   s_last(istat1),s_pts(istat1))
	       s_pas(istat1)=s_pas(istat1)+1
	       s_first(istat1)=t
	       s_pts(istat1)=1
	    endif
	    s_pts(istat1)=s_pts(istat1)+1
	    s_last(istat1)=t
	 endif

* Check if required combination occurs

	 if (mcombi.gt.0) then
	    if (t-tcombi.gt.tpass) then
	       ncombi(lcombi)=ncombi(lcombi)+1
	       lcombi=0
	    endif
	    l=0
	    do i=1,mcombi
	       if (t-s_last(combi(i)).lt.tpass) l=l+1
	    enddo
	    if (l.gt.0) then
	       lcombi=max(l,lcombi)
	       tcombi=t
	    endif
         endif
         goto 400

499      continue
	 l=index(filenm(ifile),' ')-1
	 call strf1985(date1,'%y%m%d %H:%M:%S',nint((tmin1-46066)*86400))
	 call strf1985(date2,'%y%m%d %H:%M:%S',nint((tmax1-46066)*86400))
         write (*,600) filenm(ifile)(:l),nread,date1,date2
         close (10)
      enddo

* End passes

      if (pass) then
         do istat1=1,10000
            call writepass(istat1,s_first(istat1),
     |                s_last(istat1),s_pts(istat1))
	 enddo
      endif

* Print warnings for stations with no sigma

      do istat1=1,10000
	 if (warn(istat1)) write (0,640) istat1
      enddo

* Print statistics

      close (20)
      call strf1985(date3,'%y%m%d %H:%M:%S',nint((tmin2-46066)*86400))
      call strf1985(date4,'%y%m%d %H:%M:%S',nint((tmax2-46066)*86400))
      l=index(filenm(nfile),' ')-1
      write (*,610)
     |filenm(nfile)(:l),nsel,date3,date4,ndel,nsel-ndel
600   format (/
     |'Input METRIC file : ',a/
     |'Records read      : ',i9,3x,a,' - ',a)
610   format (
     |'Output METRIC file: ',a/
     |'Records selected  : ',i9,3x,a,' - ',a/
     |'Records deleted   : ',i9/
     |'Records remaining : ',i9)
620   format (/
     |'Stat  Passes Tot.NPs Use.NPs')
630   format (/
     |'Combinations of stations',10i5)
631   format (i2,' out of',i2,':',i8)
640   format ('WARNING: no sigma defined for station',i5)

* If station table is requested

      if (stat) then
	  write (*,620)
          do istat1=1,10000
	     if (s_pnt(istat1)+s_del(istat1).gt.0)
     |		write (*,'(i4.4,3i8)') istat1,s_pas(istat1),s_pnt(istat1)
     |			+s_del(istat1),s_pnt(istat1)
          enddo
      endif

* If station combination is scanned

      if (mcombi.gt.0) then
         write (*,630) (combi(i),i=1,mcombi)
	 do i=1,mcombi
	    write (*,631) i,mcombi,ncombi(i)
	 enddo
      endif

9999  end

      subroutine writepass(i,t1,t2,n)
      integer i,n
      real*8  t1,t2,ymdhms

      if (n.le.0) return
      write (*,500)
     |	i,ymdhms((t1-46066)*86400),ymdhms((t2-46066)*86400),n
500   format (i4,2f17.3,i4)
      end
      
      subroutine readdel(file)
      character*80 file,line
      integer i,j
      real*8 a,b,sec85
      include 'gbfedit.inc'

* Open the file with delete cards and store station id and time interval

      open (10,file=file,status='old')
100   read (10,'(a)',end=199) line
      do i=41,52
         if (line(i:i).eq.' ') line(i:i)='0'
      enddo
      line(53:53)='.'
      do i=61,72
         if (line(i:i).eq.' ') line(i:i)='0'
      enddo
      line(73:73)='.'
      read (line,120) i,a,b
120   format (10x,i4,26x,2f20.7)
      if (a.eq.0) then
         a=550101000000d0	! Remove all measurements of station i
	 b=541231235959d0
      endif 
      ndel=ndel+1
      if (ndel.gt.maxdel)call fin("metricedit: too many delete cards.")
      val0(ndel)=sec85(4,a)/86400d0+46066d0
      val1(ndel)=sec85(4,b)/86400d0+46066d0
      value(ndel)=i+val0(ndel)/1d5
      pnt(ndel)=ndel
      goto 100
199   continue
      close (10)

* Now sort all cards by station id and time

      call dqsort(pnt,value,ndel)
      do i=1,ndel
         j=int(value(pnt(i)))
         delt0(i)=val0(pnt(i))
	 delt1(i)=val1(pnt(i))
	 delsta(i)=j
	 if (pnt1(j).eq.0) pnt1(j)=i
	 pnt1(j+1)=i+1
      enddo
      do i=1,10000
         pnt2(i)=pnt1(i)
      enddo	    
      end

