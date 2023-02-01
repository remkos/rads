      program giiesumm

* Program scans GEODYN II standard output files and prints out essential
* information on estimated parameters
*
* Last updates:
* 12-Jun-1997 - Added printing of Keplerian elements and sigmas
* 24-Jan-2000 - New output of station coordinates
*-----------------------------------------------------------------------
* Variable declarations

      implicit none
      integer mp,ms
      parameter (mp=500,ms=100)
      integer nbias,nbias100,nbias900,naccel,nsolrad,ndrag,npos,
     |		nsta,nsat,nstat,nper,iargc,
     |		jaccel,jsolrad,jdrag,i,ii,j,k,statyp,postyp/0/
      real*8 tdrag(mp),cdrag(mp),sdrag(mp)
      real*8 tsolrad(mp),csolrad(mp),ssolrad(mp)
      real*8 taccel(mp),caccel(mp),saccel(mp)
      real*8 cbias(mp),sbias(mp),period(ms)
      real*8 cpos(ms,18),spos(ms,18),rpos(ms,18,6)
      real*8 csta(ms,18),stavel(ms,4)
      real*8 t99,sig,t0,t1,dummy,rad,year
      parameter (year=86400d0*365.25d0)
      character*20 adrag(mp),asolrad(mp),aaccel(mp),abias(mp),
     |		apos(ms),asta(ms)
      character*132 line,line1
      integer satid900(10),satid100(10),satid(10),nc1,nc2,nc3,
     |		stanr(ms)
      character stanm(ms)*6,hpos*18
      logical outer,nonsing

* Initialise (data statements)

      data hpos /'XYZUVWAEINWMacsCSL'/
      data outer /.true./,nonsing/.false./
      data nc1/0/,nc2/0/,nc3/0/
      data nbias/0/,nbias100/0/,nbias900/0/,naccel/0/,nsolrad/0/
      data ndrag/0/,npos/0/,nsta/0/,nsat/0/,nstat/0/,nper/0/

* Initialise (remaining)

      rad=atan(1d0)/45

* Read arguments

      do i=1,iargc()
         call getarg(1,line)
	 if (line(:2).eq.'-n') then
	    nonsing=.true.
	 endif
      enddo

* Start reading INPUT DECK
* Scan for keywords:
* - station names, EPOCH, DRAG, ACCEL9, SATPAR, BIAS, ...

10    call readinput(line)
      if (line(:10).eq.'ENDGLB  11') then
         outer=.false.
      else if (line(:8).eq.'INSTRMNT') then
	 call readinput(line)
	 if (nstat.gt.ms) call fin("too many stations")
	 nstat=nstat+1
	 stanm(nstat)=line(1:6)
	 read (line(17:20),*) stanr(nstat)
*	 write (*,*) stanr(nstat),stanm(nstat)
         call readinput(line)
	 read (line(21:),*) (stavel(nstat,i),i=1,3)
         call readinput(line)
	 read (line(21:),*) stavel(nstat,4)
      else if (line(:6).eq.'EPOCH ') then
	 call readdate(line(21:40),t0)
	 call readdate(line(61:80),t0)
	 call readdate(line(41:60),t1)
      else if (line(:6).eq.'SATPAR') then
	 nsat=nsat+1
	 read (line(18:24),*) satid(nsat)
      else if (line(:4).eq.'DRAG') then
	 ndrag=ndrag+1
 	 read (line(60:72),'(d13.1)') sig
	 if (sig.le.0) then
	    ndrag=ndrag-1
	 else
	    call readdate(line(45:59),tdrag(ndrag))
	 endif
      else if (line(:6).eq.'SOLRAD') then
	 nsolrad=nsolrad+1
	 call readdate(line(45:59),tsolrad(ndrag))
 	 read (line(60:72),'(d13.1)') sig
	 if (sig.le.0) nsolrad=nsolrad-1
      else if (line(:6).eq.'ACCEL9' .and. line(11:12).eq.'99') then
	 call readdate(line(45:59),t99)
      else if (line(:6).eq.'ACCEL9') then
	 naccel=naccel+1
	 taccel(naccel)=t99
 	 read (line(60:72),'(d13.1)') sig
	 if (sig.le.0) naccel=naccel-1
      else if (line(:17).eq.'MBIAS         900') then
	 nbias900=nbias900+1
         read (line(18:24),*) satid900(nbias900)
      else if (line(:17).eq.'MBIAS         100') then
	 nbias100=nbias100+1
         read (line(18:24),*) satid100(nbias100)
      else if (line(:6).eq.'ENDALL') then
	 goto 20
      endif
      goto 10

20    continue

* Set zero time fields to t0

      do i=1,ndrag
	 if (tdrag(i).eq.0) tdrag(i)=t0
      enddo
      do i=1,nsolrad
	 if (tsolrad(i).eq.0) tsolrad(i)=t0
      enddo

* Continue scanning for parameters and sigmas till end of file

30    line1=line
      read (*,550,end=40) line
      nc1=nc1+1
      nc2=nc2+1

* Clean arc-parameter buffers when '0XPOS' or '0A   ' is found.
* Clean global-parameter buffers when '0STA' or '0LAT' is found.

      if ((line(:5).eq.'0XPOS' .or. line(:5).eq.'0A   ')
     |		.and. nc1.gt.10) then
         jaccel=0
         jdrag=0
         jsolrad=0
         nbias=0
	 nbias100=0
	 nbias900=0
	 npos=0
	 nper=0
	 nc1=0
      else if (line(:4).eq.'0STA' .and. nc2.gt.10) then
	 nsta=0
	 nc2=0
      endif

* Store principle parameters and sigmas

      if (line(:6).eq.'0STA X' .or. line(:8).eq.'0STA LAT') then
	 statyp=1
	 if (line(:8).eq.'0STA LAT') statyp=2
	 nsta=nsta+1
	 asta(nsta)=line(11:17)
	 read (line(18:),*) csta(nsta,10)
	 read (*,*) dummy,sig,csta(nsta,13)
	 read (*,*) csta(nsta,1),dummy,csta(nsta,4)
	 if (sig.eq.0) csta(nsta,4)=-1
	 nc2=0
      else if (line(:6).eq.'0STA Y' .or. line(:8).eq.'0STA LON') then
	 read (line(18:),*) csta(nsta,11)
	 read (*,*) dummy,sig,csta(nsta,14)
	 read (*,*) csta(nsta,2),dummy,csta(nsta,5)
	 nc2=0
      else if (line(:6).eq.'0STA Z' .or. line(:8).eq.'0STA HGH') then
	 read (line(18:),*) csta(nsta,12)
	 read (*,*) dummy,sig,csta(nsta,15)
	 read (*,*) csta(nsta,3),dummy,csta(nsta,6)
	 nc2=0
      else if (line(:40).eq.'     EARTH FIXED RECTANGULAR COORDINATES')
     |then
         nsta=0
32       read (*,550,end=40) line
         if (line(:9).eq.'  APRIORI') then
	    nsta=nsta+1
	    read (line(102:),*) (csta(nsta,i),i=16,18)
	 else if (line(:9).eq.' ADJUSTED') then
	    read (line(102:),*) (csta(nsta,i),i=7,9)
	 else if (line(:1).eq.'0') then
	    goto 34
	 endif
	 goto 32
34       nc2=0
      else if (line(:5).eq.'0XPOS') then
	 npos=npos+1
	 apos(npos)=line(11:17)
	 call skip(1)
	 read (*,*) cpos(npos,1),sig,spos(npos,1)
	 nc1=0
	 nc3=0
	 postyp=1
      else if (line(:5).eq.'0YPOS') then
	 call skip(1)
	 read (*,*) cpos(npos,2),sig,spos(npos,2)
	 nc1=0
      else if (line(:5).eq.'0ZPOS') then
	 call skip(1)
	 read (*,*) cpos(npos,3),sig,spos(npos,3)
	 nc1=0
      else if (line(:5).eq.'0XVEL') then
	 call skip(1)
	 read (*,*) cpos(npos,4),sig,spos(npos,4)
	 nc1=0
      else if (line(:5).eq.'0YVEL') then
	 call skip(1)
	 read (*,*) cpos(npos,5),sig,spos(npos,5)
	 nc1=0
      else if (line(:5).eq.'0ZVEL') then
	 call skip(1)
	 read (*,*) cpos(npos,6),sig,spos(npos,6)
	 nc1=0
      else if (line(:5).eq.'0A   ') then
	 npos=npos+1
	 apos(npos)=line(11:17)
	 call skip(1)
	 read (*,*) cpos(npos,7),sig,spos(npos,7)
	 nc1=0
	 nc3=0
	 postyp=2
      else if (line(:5).eq.'0ECC ') then
	 call skip(1)
	 read (*,*) cpos(npos,8),sig,spos(npos,8)
	 nc1=0
      else if (line(:5).eq.'0INCL') then
	 call skip(1)
	 read (*,*) cpos(npos,9),sig,spos(npos,9)
	 nc1=0
      else if (line(:5).eq.'0R A ') then
	 call skip(1)
	 read (*,*) cpos(npos,10),sig,spos(npos,10)
	 nc1=0
      else if (line(:5).eq.'0ARG ') then
	 call skip(1)
	 read (*,*) cpos(npos,11),sig,spos(npos,11)
	 nc1=0
      else if (line(:5).eq.'0M AN') then
	 call skip(1)
	 read (*,*) cpos(npos,12),sig,spos(npos,12)
	 nc1=0
      else if (index(line,'  RECTANGULAR').gt.0
     |  .and. postyp.ne.1) then
         call skip(9)
	 read (*,550) line
	 read (line(27:),*) (cpos(npos,j),j=1,6)
	 call skip(5)
	 read (*,550) line
	 read (line(27:),*) (spos(npos,j),j=1,6)
         nc3=0
      else if (index(line,'  KEPLERIAN').gt.0
     |  .and. postyp.ne.2) then
         call skip(9)
	 read (*,550) line
	 read (line(27:),*) (cpos(npos,j),j=7,12)
	 call skip(5)
	 read (*,550) line
	 read (line(27:),*) (spos(npos,j),j=7,12)
         nc3=0
      else if
     |  (index(line,'  NON SINGULAR').gt.0
     |  .and. postyp.ne.3) then
         call skip(9)
	 read (*,550) line
	 read (line(27:),*) (cpos(npos,j),j=13,18)
	 call skip(5)
	 read (*,550) line
	 read (line(27:),*) (spos(npos,j),j=13,18)
         nc3=0
      else if (index(line,'ADJUSTED CARTESIAN').gt.0
     |	.and. (line1.eq.' ' .or. line1(:5).eq.' DRAG')
     |  .and. nc3.eq.0 .and. .not.outer) then
         call skip(1)
	 read (*,*) dummy,(rpos(npos,1,j),j=2,6)
	 read (*,*) dummy,(rpos(npos,2,j),j=3,6)
	 read (*,*) dummy,(rpos(npos,3,j),j=4,6)
	 read (*,*) dummy,(rpos(npos,4,j),j=5,6)
	 read (*,*) dummy,(rpos(npos,5,j),j=6,6)
	 do j=1,6
	    do k=1,j-1
	       rpos(npos,j,k)=rpos(npos,k,j)
	    enddo
	    rpos(npos,j,k)=1
	 enddo
	 nc3=1
      else if (index(line,'ADJUSTED KEPLERIAN').gt.0
     |	.and. (line1.eq.' ' .or. line1(:5).eq.' DRAG')
     |  .and. nc3.eq.0 .and. .not.outer) then
         call skip(1)
	 read (*,*) dummy,(rpos(npos,7,j),j=2,6)
	 read (*,*) dummy,(rpos(npos,8,j),j=3,6)
	 read (*,*) dummy,(rpos(npos,9,j),j=4,6)
	 read (*,*) dummy,(rpos(npos,10,j),j=5,6)
	 read (*,*) dummy,(rpos(npos,11,j),j=6,6)
	 do j=1,6
	    do k=1,j-1
	       rpos(npos,j+6,k)=rpos(npos,k+6,j)
	    enddo
	    rpos(npos,j+6,k)=1
	 enddo
	 nc3=1
      else if (index(line,'ADJUSTED NON SINGULAR').gt.0
     |	.and. (line1.eq.' ' .or. line1(:5).eq.' DRAG')
     |  .and. nc3.eq.0 .and. .not.outer) then
         call skip(1)
	 read (*,*) dummy,(rpos(npos,13,j),j=2,6)
	 read (*,*) dummy,(rpos(npos,14,j),j=3,6)
	 read (*,*) dummy,(rpos(npos,15,j),j=4,6)
	 read (*,*) dummy,(rpos(npos,16,j),j=5,6)
	 read (*,*) dummy,(rpos(npos,17,j),j=6,6)
	 do j=1,6
	    do k=1,j-1
	       rpos(npos,j+12,k)=rpos(npos,k+12,j)
	    enddo
	    rpos(npos,j+12,k)=1
	 enddo
	 nc3=1
      else if (line(:3).eq.'0GA') then
	 jaccel=jaccel+1
	 aaccel(jaccel)=line(8:17)
	 if (aaccel(jaccel)(:2).eq.'99') then
	    if (aaccel(jaccel-1)(:2).ne.'23') then
	       write (0,*) 'something fishy with ACCEL9 99 card'
	    endif
	    aaccel(jaccel)(:2)='33'
	 endif
	 call skip(1)
	 read (*,*) caccel(jaccel),sig,saccel(jaccel)
      else if (line(:3).eq.'0CD') then
	 jdrag=jdrag+1
	 adrag(jdrag)=line(11:17)
	 call skip(1)
	 read (*,*) cdrag(jdrag),sig,sdrag(jdrag)
      else if (line(:3).eq.'0CR') then
	 jsolrad=jsolrad+1
	 asolrad(jsolrad)=line(11:17)
	 call skip(1)
	 read (*,*) csolrad(jsolrad),sig,ssolrad(jdrag)
      else if (line(:2).eq.'0B') then
	 nbias=nbias+1
	 if (line(10:16).ne.' ') then
	 else if (line(3:5).eq.'900') then
	    nbias900=nbias900+1
	    write (line(10:16),'(i7)') satid900(nbias900)
	 else if (line(3:5).eq.'100') then
	    nbias100=nbias100+1
	    write (line(10:16),'(i7)') satid100(nbias100)
	 else 
	    read (line(7:9),*) i
	    write (line(10:16),'(i7)') satid(i+1)
*	    write (*,*) 'satid',i
	 endif
	 abias(nbias)=line(3:16)
	 call skip(1)
	 read (*,550) line
	 read (line,*) cbias(nbias),sig,sbias(nbias)
      else if (line(:9).eq.' PERIOD =') then
	 nper=nper+1
         read (line(10:29),*) period(nper)
	 period(nper)=period(nper)*60
      endif
      goto 30

40    continue

* Print out results

      if(jaccel.ne.naccel) write (0,*) 'something wrong with ACCEL9',
     |	jaccel,naccel
      if(jdrag.ne.ndrag) write (0,*) 'something wrong with DRAG',
     |	jdrag,ndrag
      if(jsolrad.ne.nsolrad)write (0,*)'something wrong with SOLRAD',
     |	jsolrad,nsolrad

550   format (a)
600   format(a18,t20,f18.3,d18.12,d18.12)
605   format(a18,t20,2f18.3)
610   format(a18,t20,3f18.6,6f12.6)
615   format(a18,t20,3f18.9,f12.3)
620   format(a18,t20,3f18.9,6f12.6)
630   format(a18,t20,f18.6,f18.12,f18.9)
640   format(a18,t20,f18.6,2f18.12)
650   format(a18,t20,2f18.12,f18.9)
660   format(a18,t20,6f8.4)

      write (*,605) 'EPOCH             ',t0,t1
      do i=1,nsta
	 if (csta(i,4).ge.0) then
	 write(*,610)'STA        '//asta(i),(csta(i,j),j=1,9)
	 write(*,610)'STAAPR     '//asta(i),(csta(i,j),j=10,18)
	 read (asta(i),*) j
	 do k=1,nstat
	    if (j.eq.stanr(k)) goto 101
	 enddo
	 call fin('Could not find station nr')
101	 write(*,615)'STAVEL     '//asta(i),(stavel(k,j)*year,j=1,3),
     |		stavel(k,4)
         endif
      enddo

      do i=1,npos
	 write(*,610)'POS        '//apos(i),(cpos(i,j),j=1,3),
     |	     (spos(i,j),j=1,3)
	 write(*,620)'VEL        '//apos(i),(cpos(i,j),j=4,6),
     |	     (spos(i,j),j=4,6)
	 write(*,630)'AEI        '//apos(i),(cpos(i,j),j=7,9)
	 write(*,630)'AEISIG     '//apos(i),(spos(i,j),j=7,9)
	 write(*,620)'NWM        '//apos(i),(cpos(i,j),j=10,12)
	 write(*,620)'NWMSIG     '//apos(i),(spos(i,j),j=10,12)
	 if (nonsing) then
	    write(*,640)'acs        '//apos(i),(cpos(i,j),j=13,15)
	    write(*,640)'acsSIG     '//apos(i),(spos(i,j),j=13,15)
	    write(*,650)'CSL        '//apos(i),(cpos(i,j),j=16,18)
	    write(*,650)'CSLSIG     '//apos(i),(spos(i,j),j=16,18)
	    ii=18
	 else
	    ii=12
	 endif
	 if (outer) ii=0
	 do j=1,ii
	    write(*,660)'CORREL '//hpos(j:j)//'   '//apos(i),
     |	     (rpos(i,j,k),k=1,6)
	 enddo
	 write(*,620)'PERIOD     '//apos(i),period(i)
      enddo

      do i=1,ndrag
         write(*,600)'DRAG       '//adrag(i),tdrag(i),cdrag(i),
     |		sdrag(i)
      enddo
      do i=1,nsolrad
         write(*,600)'SOLRAD     '//asolrad(i),tsolrad(i),csolrad(i),
     |		ssolrad(i)
      enddo
      do i=1,naccel
         write(*,600)'ACCEL9  '//aaccel(i),taccel(i),caccel(i),
     |		saccel(i)
      enddo
      do i=1,nbias
	 if (cbias(i).ne.0d0) then
	    do j=1,nstat
	       if (abias(i)(8:13).eq.stanm(j)) then
		  write (abias(i)(8:14),'(i4,3x)') stanr(j)
	       endif
	    enddo
	    write(*,600)'BIAS   '//
     |		abias(i)(:3)//' '//abias(i)(8:14),
     |		t0,cbias(i),sbias(i)
	 endif
      enddo
      end

      subroutine readdate(text,date)
      character*14 text
      real*8 date
      integer i
      do i=1,14
	 if (text(i:i).eq.' ') text(i:i)='0'
      enddo
      text(13:13)='.'
      read (text,'(f14.1)',err=9) date
      goto 10

   9  date=0.0D0

  10  continue

      end

      subroutine readinput(line)
      character*(*) line
      character*86 text

10    read (*,'(a)',end=999) text
      if (index(text,'EXECUTION TERMINATED').gt.0) stop
      if (text(81:86).eq.' ') goto 10
      line=text(2:)
999   end

      subroutine skip(n)
      integer i,n
      do i=1,n
         read (*,*)
      enddo
      end
