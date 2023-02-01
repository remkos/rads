      program g2todr

* Convert GEODYN II Trajectory file (G2T) to Orbital Data Records
* (ODR). The can also update the arclist and archive the ODR and arclist
* in the proper directory and send it to NOAA.
*
* To list the content of the resulting ODR, use lodr.
*-
* $Log: g2todr.f,v $
* Revision 1.4  2007/02/21 18:26:30  rads
* - Small changes to work better with gfortran (gcc 4.2)
*
* Revision 1.3  2004/11/22 14:34:09  remko
* - Allow arc numbers higher than 999
*
* 18-Oct-2002 - Use iargc for compatibility with Intel f95 compiler
* 13-Sep-2002 - Make high-precision ODR format default, new +x flag
*  2-May-2002 - Handles little endian machines properly
*  8-Jan-1999 - Added EPOCH card with ELEM option (-te)
* 21-Oct-1996 - Supports new high-precision ODR format
* 27-Jan-1992 - First version
************************************************************************
      implicit none

      integer      maxmod,maxrec
      parameter    (maxmod=25,maxrec=40000)
      real*8	   buffer(2048)
      character*80 deck(250),filenm,arg,modelnm,line,parname,
     |		   resname,arclistname,shellname,odrname
      character    spec*4,satel*8,model(maxmod)*15,
     |		   date0*15,date1*15,date2*15,dirnames*160
      integer*4    irem(maxmod),m_ndeg(maxmod),pole/0/,rem,rep,arc/0/
      logical      local/.false./,noaa/.false./,plain/.false./,
     |		   batch/.false./,greenw/.false./,elems/.false./,
     |		   buff/.false./,exist/.false./,newformat/.true./,
     |		   ltlend
      real*8       mjdarc0,mjdint,geo(3),xyz(6)
      real*8	   ecf(3),m_ae(maxmod),tskip/-1/
      integer	   isat,nsat,maxsat
      parameter    (maxsat=50)
      integer	   nrec,isatel,satid(maxsat),ios,te/0/
      equivalence  (buffer(49),deck)
      real*8       pi,rad,masrad,scale/1d6/,sec85
      integer*4    itime,lat,lon,hgt,odr4(4),head1,head2(4)
      integer	   i,j,krec,idate,imin,ihour,idat,ibuf,ntb,
     |		   itime0/999999999/,itime1/-999999999/,
     |		   nsta,npnt,ntimbf,ioff,ndeg,nwdsat,mdate,mjd85,iarg,na,
     |		   leapsec,mjd,imodel,nrms,iarc/0/,iskip,nmodel,irep
      real*8	   grhran,xpole,ypole,time,flat,gm,ae,dtime,skip,dlon,
     |		   dlat,dh,rms,r

      data	   rem/-1/,rep/-1/,resname/' '/,arclistname/' '/,
     |		   shellname/' '/,odrname/' '/,satel/'SATELLIT'/
      equivalence  (odr4(1),itime),(odr4(2),lat),(odr4(3),lon),
     |		   (odr4(4),hgt)

      nsat=0
      skip=0
      nrec=0
      mjd85=mdate(2,850101)
      filenm='fort.30'
      pi=4*atan(1d0)
      rad=pi/180
      masrad=rad/3600d3

* Scan the command line

      iarg=0
   10 call getnextarg(iarg,arg)
      if (arg(:2).eq.'-o') then
	 plain=.true.
	 batch=.true.
	 odrname=arg(3:)
	 if (odrname.eq.' ') call getnextarg(iarg,odrname)
	 if (odrname.eq.' ') odrname='ODR.tmp'
      else if (arg(:2).eq.'n=') then
	 read (arg(3:),*) nrec
      else if (arg(:2).eq.'t=') then
	 read (arg(3:),*) tskip
      else if (arg(:3).eq.'-te') then
	 elems=.true.
	 r=0
	 te=0
	 read (arg(4:),*,iostat=ios) r
	 if (r.ne.0) te=nint(sec85(0,r))
      else if (arg(:3).eq.'-tg') then
	 greenw=.true.
      else if (arg(:2).eq.'-a') then
	 local=.true.
	 arclistname='arclist.tmp'
	 noaa=.true.
      else if (arg(:2).eq.'-l') then
	 local=.true.
      else if (arg(:2).eq.'-u') then
	 arclistname='arclist.tmp'
      else if (arg(:2).eq.'-n') then
	 noaa=.true.
      else if (arg(:2).eq.'-S') then
	 shellname=arg(3:)
	 if (shellname.eq.' ') shellname='shell.tmp'
	 batch=.true.
      else if (arg(:2).eq.'-R') then
	 resname=arg(3:)
      else if (arg(:2).eq.'-A') then
	 arclistname=arg(3:)
	 if (arclistname.eq.' ') arclistname='arclist.tmp'
      else if (arg(:2).eq.'-b') then
	 batch=.true.
      else if (arg(:2).eq.'-p') then
	 pole=1
      else if (arg(:2).eq.'+p') then
	 pole=-1
      else if (arg(:2).eq.'-r') then
      else if (arg(:2).eq.'-x') then
         newformat=.true.
      else if (arg(:2).eq.'+x') then
         newformat=.false.
      else if (arg(:2).eq.'-B') then
	 buff=.true.
      else if (arg(:2).eq.'-h') then
	 call system('cat /user/altim/src/gdyn/g2todr.man')
	 goto 9999
      else if (arg(:4).eq.'sat=') then
         satel=arg(5:)
      else if (arg(:6).eq.'model=') then
	 modelnm=arg(7:)
      else if (arg(:4).eq.'arc=') then
         read (arg(5:),'(z1,i2)') i,j
	 arc=i*100+j
      else if (arg(:4).eq.'rem=') then
	 read (arg(5:),*) rem
      else if (arg(:4).eq.'rep=') then
	 read (arg(5:),*) rep
      else if (arg(:1).eq.' ') then
	 goto 15
      else
	 filenm=arg
	 modelnm=filenm
      endif
      goto 10

15    continue
      if (local .or. noaa .or. plain) then
	 if (newformat) then
	    scale=1d7
	    spec='xODR'
	 else
	    scale=1d6
	    spec='@ODR'
	 endif
      else
	 spec=' '
      endif

* Open G2T file

      open (30,status='old',file=filenm,form='unformatted',err=1350)

* Read G2T header buffer

      rewind (30)
      read (30) buffer
      if (buffer(1).ne.-9d9) write (0,*) 'buffer(1)=',buffer(1)
      na=nint(buffer(2))
      if (buffer(7).ne.1) goto 1340
      nwdsat=nint(buffer(8))
      ntimbf=nint(buffer(10))
      dtime=buffer(19)
      ioff=0
      do 20 i=202,207
   20    if (buffer(i).gt.0) ioff=ioff+1
      if (buffer(208).le.0 .or. buffer(209).le.0 .or.
     |    buffer(210).le.0) goto 1320
      isat=nint(buffer(301))
      if (isat.gt.maxsat) call fin("g2todr: too many satellites")

* Read Alphanumeric header(s)

      do j=1,na
	 read (30,end=1330) buffer
	 do i=1,250
	    line=deck(i)
	    if (line(:6).eq.'SATPAR') then
	       nsat=nsat+1
	       if (nsat.gt.maxsat)call fin("g2tpoe: too many satellites")
	       read (line(18:24),'(bz,i7.7)') satid(nsat)
	    else if (line(:5).eq.'EARTH') then
	       read (line,600) ndeg,gm,ae,flat
  600 format (14x,i3,7x,d20.8,d15.3,d13.1)
	    endif
	 enddo
         if (buffer(1).ne.-8d9) backspace (30)
      enddo

* Open ODR-file if requested

      if (spec.ne.' ')
     |	open (20,file=odrname,status='new',form='unformatted',
     |	recl=16,access='direct')

      if (pole.eq.0 .and. gm.ne.398600.436d9) pole=1

* Process Data buffers

      time=0
      krec=0
  100 read (30,end=1330) buffer
      if (buffer(1).eq.9d9) goto 200

* Store date and start-time

      idate=int(buffer(2)/1d6)
      time=buffer(2)-idate*1d6
      ihour=int(time/1d4)
      time=time-ihour*1d4
      imin=int(time/1d2)
      time=time-imin*1d2
      time=ihour*3600+imin*60+time
      ntb=nint(buffer(5))
      ibuf=2*ntimbf-1
      do idat=1,ntb
	 do i=1,6
	    xyz(i)=buffer(ibuf+ioff+i)
	 enddo
	 do i=1,3
	    geo(i)=buffer(ibuf+ioff+6+i)
	 enddo
	 grhran=buffer(5+ntimbf+idat)
	 xpole=buffer(ibuf+ioff+16)*masrad
	 ypole=buffer(ibuf+ioff+17)*masrad
	 if (pole.eq.1) then
	    xpole=buffer(ibuf+ioff+16)*masrad
	    ypole=buffer(ibuf+ioff+17)*masrad
	 else
	    xpole=0
	    ypole=0
	 endif
	 dlat=geo(1)*rad
	 dlon=geo(2)*rad
	 dh=-ae/flat*sin(2*dlat)*(xpole*cos(dlon)-ypole*sin(dlon))
	 call rotate(3,grhran,xyz,ecf)
	 call rotate(1,-ypole,ecf,ecf)
	 call rotate(2,-xpole,ecf,ecf)
* write (*,*) xyz,ecf,xpole,ypole
	 call xyzgeo(ecf,r,geo(1),geo(2),geo(3))
	 geo(1)=geo(1)/rad
	 geo(2)=geo(2)/rad

* For new ODR format: -180 < LON <= 180
* Otherwise:          0 <= LON < 360

	 if (spec.eq.'xODR') then
            if (geo(2).gt.180) geo(2)=geo(2)-360
         else
	    if (geo(2).lt.0) geo(2)=geo(2)+360
         endif
	 krec=krec+1
	 if (krec.gt.maxrec)call fin("g2todr: too many records")
	 mjd=mdate(2,idate)-mjd85
	 itime=mjd*86400+nint(time)
	 itime0=min(itime0,itime)
	 itime1=max(itime1,itime)
	 if (spec.ne.' ') then
	    lat=nint(geo(1)*scale)
	    lon=nint(geo(2)*scale)
	    hgt=nint(geo(3)*1d3)
	    if (ltlend()) call i4swap(4,odr4)
	    write (20,rec=krec+2) odr4
	 else
	    call strf1985(date0,'%y%m%d%H%M%S',itime)
	    if (elems) then
     	       if (te.eq.0 .or. abs(te-itime).le.10)
     |		write (*,640) date0,xyz
	    else if (greenw) then
	       write (*,630) date0,geo,grhran*rad
	    else
	       write (*,630) date0,geo
	       if (buff) write (*,650) (buffer(ibuf+ioff+i),i=1,20),dh
	    endif
	 endif
	 ibuf=ibuf+nwdsat
         time=time+dtime
      enddo
      goto 100

200   if (spec.eq.' ') goto 9999

* Check if there was a leapsecond in the data period
* If so, determine when the leapsecond occurred, and correct all
* timetags after the event.
* This is because of the block structure of the G2T file.

      if (mod(itime0,10).eq.0 .and. mod(itime1,10).eq.9) then
	 mjd=itime1/86400+mjd85
	 mjd=mdate(2,mdate(1,mjd)/100*100+01)
	 leapsec=(mjd-mjd85)*86400
	 do i=3,krec+2
	    read (20,rec=i) odr4
	    if (ltlend()) call i4swap(1,itime)
	    if (itime.gt.leapsec .and. mod(itime,10).eq.0) then
	       itime=itime-1
	       if (ltlend()) call i4swap(1,itime)
	       write (20,rec=i) odr4
	    endif
	 enddo

* Another strange thing occurs when the leapsecond IS the very first
* second. It sets the first block one second too far.
* Thus we correct all records with a minute+1second to the entire minute

      else if (mod(itime0,10).eq.1 .and. mod(itime1,10).eq.0) then
         do i=3,krec+2
	    read (20,rec=i) odr4
	    if (ltlend()) call i4swap(1,itime)
	    if (mod(itime,10).eq.1) then
	       itime=itime-1
	       if (ltlend()) call i4swap(1,itime)
	       write (20,rec=i) odr4
	    endif
	 enddo

* Now if all is not well, we give an error.

      else if (mod(itime0,10).ne.0 .or. mod(itime1,10).ne.0) then
         call fin("g2todr: time tags are not on the entire 10 seconds")
      endif

* Get the parameter file

      parname='/user/altim'
      call checkenv('ALTIM',parname,i)
      parname(i+1:)='/src/gdyn/g2todr.par'
      open (9,file=parname,status='old')
      rewind (9)
  205 read (9,*,iostat=i) arg
      if (i.ne.0) then
         write (*,1361) satid(isat)
 1361 format ('g2todr: end of parameter file:',
     |' satid',i8,' not found. Continuing.')
	 irep=rep
         goto 206
      endif
      read (9,*) isatel
      read (9,*) dirnames
      read (9,*) mjdarc0,mjdint,skip
      read (9,*) irep
      read (9,*) nmodel
      do i=1,nmodel
	 read (9,*) m_ndeg(i),m_ae(i),model(i),irem(i)
	 if (rem.ne.-1) irem(i)=rem
      enddo
      read (9,*) nsta
      read (9,*) npnt
      read (9,*) nrms
      if (isatel.ne.satid(isat)) goto 205
      satel=arg(:8)
      iarc=nint((itime0/86400d0+mjd85-mjdarc0)/mjdint)
      if (rep.ne.-1) then
	 irep=rep
      else if (irep.ne.-1) then
      else if (iarc.lt.76) then
	 irep=3000
      else if (iarc.eq.76) then
	 irep=0
      else if (iarc.lt.256) then
	 irep=35000
      else if (iarc.eq.256) then
	 irep=0
      else if (iarc.lt.286) then
	 irep=3000
      else if (iarc.eq.286) then
	 irep=0
      else
	 irep=168000
      endif
206   if (arc.ne.0) iarc=arc

* Test for model

      do imodel=1,nmodel
	 if (m_ndeg(imodel).eq.ndeg .and. m_ae(imodel).eq.ae) goto 210
      enddo
      imodel=1
  210 continue
      if (odrname.eq.' ')write(odrname,660)'ODR.',iarc/100,mod(iarc,100)

* Can you find a residual file ?

      if (.not.batch) then
      if (resname.eq.' ')write(resname,660)'res.',iarc/100,mod(iarc,100)
      inquire (exist=exist,file=resname)
      if (.not.exist) resname="lres.out"
      inquire (exist=exist,file=resname)
      if (exist) then
	 open (50,file=resname,status='old')
	 rewind (50)
  220	 read (50,550,end=230) arg
	 goto 220
  230    read (arg(37:41),*) nsta
	 read (arg(56:62),*) npnt
	 read (arg(70:80),*) rms
	 nrms=nint(rms)
      endif

* Test what has been inputted

      write (*,552) 'Gravity model',model(imodel)
      read (0,550) arg
      if (arg.ne.' ') model(imodel)=arg(:15)
      write (*,552) 'Satellite',satel
      read (0,550,end=9999) arg
      if (arg.ne.' ') satel=arg(:8)
      write (*,553) 'Repeat ID',irep
      read (0,550,end=9999) arg
      if (arg.ne.' ') read (arg,620) irep
      write (*,553) 'Arc number',iarc
      read (0,550,end=9999) arg
      if (arg.ne.' ') read (arg,620) iarc
      write (*,553) 'Remark ID',irem(imodel)
      read (0,550,end=9999) arg
      if (arg.ne.' ') read (arg,620) irem(imodel)
      write (*,552) 'ODR filename',odrname
      read (0,550,end=9999) arg
      if (arg.ne.' ') odrname=arg
      if (arclistname.ne.' ') then
      write (*,553) 'Number of stations',nsta
      read (0,550,end=9999) arg
      if (arg.ne.' ') read (arg,620) nsta
      write (*,553) 'Number of SLR data',npnt
      read (0,550,end=9999) arg
      if (arg.ne.' ') read (arg,620) npnt
      write (*,553) 'RMS of fit in cm',nrms
      read (0,550,end=9999) arg
      if (arg.ne.' ') read (arg,620) nrms
      endif
      endif

* Search for begin of precise arc

      skip=itime0+skip*86400
      call search(dble(itime0),60d0,skip,krec)
      iskip=nint(skip)
      if (tskip.ne.-1) iskip=nint(tskip)
      if (nrec.ne.0)  krec=nrec

* Write headers

      head1   =iskip
      head2(1)=irep
      head2(2)=iarc
      head2(3)=krec
      head2(4)=irem(imodel)
      if (ltlend()) then
         call i4swap(1,head1)
         call i4swap(4,head2)
      endif
      write (20,rec=1) spec,satel,head1
      write (20,rec=2) head2
      close (20)

* Update arclist

      if (arclistname.ne.' ') then
	 open (21,file=arclistname,status='new')
	 call strf1985(date0,'%y%m%d %H:%M',itime0)
	 call strf1985(date1,'%y%m%d %H:%M',itime1)
	 call strf1985(date2,'%y%m%d %H:%M:%S',iskip )
	 write (21,670) iarc/100,mod(iarc,100),date0(:12),date1(:12),
     |		nsta,npnt,nrms,irep/1e3,irem(imodel),date2
	 close (21)
      endif

* Write the environment variables

      if (shellname.ne.' ') then
	 open (21,file=shellname,status='new')
	 write (21,550) "setenv model "//model(imodel)
	 write (21,550) "setenv odrname "//odrname
	 write (21,550) "setenv arclistname "//arclistname
	 write (21,550) "setenv resname "//resname
	 write (21,660) "setenv arc ",iarc/100,mod(iarc,100)
	 write (21,550) dirnames
	 close (21)
      endif

* Move to right directory

      if (.not.plain)
     |call system('/home/dut/vlrusch/bin/odr2dutrex')

* Send to local

      if (local)
     |call system('/home/dut/vlrusch/bin/odr2ssrt')

* Send to NOAA

      if (noaa)
     |call system('/home/dut/vlrusch/bin/odr2noaa')

      goto 9999

* Formats

  550 format (a)
  552 format (a,' [',a,'] -> ',$)
  553 format (a,' [',i6,'] -> ',$)
* 600 format ('G2T: ',f15.1,-9pf15.4,0p,2f15.4,2i3)
* 610 format (a/a/a)
  620 format (i60)
  630 format (a15,2f13.7,f12.3,2f14.9)
  640 format ('EPOCH               ',a12,
     |'.00000  YYMMDDHHMMSA.00000  YYMMDDHHMMSC.00000'/
     |'ELEMS110   300  ',3f20.6/'ELEMS2',10x,3f20.9)
  650 format (5f15.4)
  660 format (a,z1.1,i2.2)
  670 format (z1.1,i2.2,2x,a,' - ',a,2i5,i6,f9.3,i4,2x,a)

* Errors

 1320 write (*,550)
     |'g2todr: latitude, longitude or height not available'
      goto 9999

 1330 write (*,550) 'g2todr: premature end of file'
      goto 9999

 1340 write (*,550) 'g2todr: multiple satellites not allowed'
      goto 9999

 1350 write (*,550) 'g2todr: can''t open G2T file'
      goto 9999

 9999 end

      subroutine getnextarg(iarg,arg)
      implicit none
      character*(*) arg
      integer iarg, iargc
      iarg=iarg+1
      if (iarg.gt.iargc()) then
        arg = ' '
      else
        call getarg(iarg,arg)
      endif
      end

      subroutine search(time0,step,skip,krec)
      implicit none
      real*8 time0,step,skip
      integer krec,irec
      integer itime0,itime1,itime2,lat0,lat1,lat2
      real*8 dlat,x
      logical ltlend

      irec=min(nint((skip-time0)/step)+1,krec-1)
   10 read (20,rec=2+irec-1) itime0,lat0
      read (20,rec=2+irec  ) itime1,lat1
      read (20,rec=2+irec+1) itime2,lat2
      if (ltlend()) then
         call i4swap(1,itime0)
         call i4swap(1,lat0)
         call i4swap(1,itime1)
         call i4swap(1,lat1)
         call i4swap(1,itime2)
         call i4swap(1,lat2)
      endif
      if (lat1.lt.lat0 .and. lat1.lt.lat2 .or. irec.eq.krec-1) then
         call inter(dble(lat0),dble(lat1),dble(lat2),x,dlat)
c write (*,*) lat0,lat1,lat2,int(dlat),x
         skip=itime1+x*step
      else if (lat2.lt.lat0) then
         irec=irec+1
         goto 10
      else
         irec=irec-1
         goto 10
      endif
      end

      subroutine inter(x0,x1,x2,t,xt)
      implicit none
      real*8 x0,x1,x2,t,xt,b,c
      b=(x2-x0)/2
      c=x2-b-x1
      t=-b/(2*c)
      xt=x1+b*t+c*t*t
      end
