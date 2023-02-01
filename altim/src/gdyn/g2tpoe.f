      program g2tlist
************************************************************************
*     This programme converts data records from g2t to poe format      *
*     Don't forget to link the library " ~vlrusch/libsom/rssubs.a "    *
************************************************************************
      implicit none
      character*80 deck(250),filenm/' '/,arg,line
      character*6  program,obstyp/' '/,obslev/' '/
      character*150 prcrec /' '/
      character    satel*4,ctype*3,orbtyp*1
      real*8       buffer(2048),geo(3),xyz(6),ecf(6),pi,degrad,
     .	dtime,time,time0,grhran,xpole,ypole,et,sec,r,z/1d0/,
     .  secday/86400d0/,et0,et1,et2000,xpole0/45d0/,ypole0/286d0/,
     .  masrad,pitch,roll,yaw,u,elem(6),gm,xmjd,xmjd0,
     .  xmjdp0/0d0/,zdot/0d0/
      equivalence  (buffer(49),deck)
      integer	   isat,nsat/0/,maxsat
      parameter    (maxsat=5)
      integer      type,idate,mdate,mjd,iyear,imonth,iday,ihour,imin,
     .	na,nwdsat,ntimbf,ioff,ibuf,i,ntb,idat,mjd0,jdyear,
     .  satid(maxsat),j,ascarc,modid,relid,qualit,quali/2/,iarg,iargc
      logical	   allarc

      integer      atc1(22)
      real*8       atc2(4)

      allarc = .true.
      pi=4*atan(1d0)
      degrad=pi/180		! convert deg -> rad by multiplying with degrad
      masrad=degrad/3600d3	! convert mas -> rad by multiplying with masrad
      jdyear=0
      et2000=mdate(2,000101)+0.5d0-30000d0 ! ET days of 1.5 Jan 2000

* Set dummy attitude control record zero

      do 11 i=1,22
         atc1(i)=0
11    continue
      do 12 i=1,4
         atc2(i)=0
12    continue

* Which program ?

      call getarg(0,arg)
      if (index(arg,'g2tpoe').ne.0) then
	 program='g2tpoe'
	 type=1
      else if (index(arg,'g2tpon').ne.0) then
	 program='g2tpon'
	 type=2
      else if (index(arg,'g2tdut').ne.0) then
	 program='g2tdut'
	 type=3
      else if (index(arg,'g2tprc').ne.0) then
	 program='g2tprc'
	 ctype='PRC'
	 orbtyp='P'
	 type=4
      else if (index(arg,'g2tprl').ne.0) then
	 program='g2tprl'
	 ctype='PRL'
	 orbtyp='V'
	 type=4
      else
	 goto 1360
      endif

* Read arguments

      filenm=' '
      do iarg=1,iargc()
	 call getarg(iarg,arg)
	 if (arg(1:2).eq.'-a') then
	    allarc=.true.
	 else
	    filenm=arg
	 endif
      enddo

* Open G2T file

      if (filenm.eq.' ') filenm='fort.30'
      open (30,status='old',file=filenm,form='unformatted',err=1350)

* Read G2T header buffer

10    rewind (30)
      read (30,err=1300) buffer
      if (buffer(1).ne.-9d9) goto 1310
      na=nint(buffer(2))
      if (buffer(7).ne.1) goto 1340
      nwdsat=nint(buffer(8))
      ntimbf=nint(buffer(10))
      et0=(buffer(15)+buffer(16))/secday
      et1=(buffer(17)+buffer(18))/secday
      dtime=buffer(19)
      ioff=0
      do 20 i=202,207
   20    if (buffer(i).gt.0) ioff=ioff+1
      if (buffer(208).le.0 .or. buffer(209).le.0 .or.
     .    buffer(210).le.0) goto 1320
      isat=nint(buffer(301))
      if (isat.gt.maxsat) call fin("g2tpoe: too many satellites")

* Read Alphanumeric header(s)

      do j=1,na
         read (30,err=1300,end=1330) buffer
         do i=1,250
	    line=deck(i)
	    if (line(1:6).eq.'SATPAR') then
	       nsat=nsat+1
	       if (nsat.gt.maxsat) call fin("g2tpoe: too many satellites")
	       read (line(18:24),'(bz,i7)') satid(nsat)
	    else if (line(1:5).eq.'ACCEL')then
	       qualit=qualit+1
	    else if (line(1:5).eq.'SIGMA' .and. line(16:17).eq.'51') then
	       obstyp(1:2)='LA'
	       obslev(1:2)='QL'
	    else if (line(1:5).eq.'OBSCO' .and. line(16:17).eq.'51') then
	       obstyp(1:2)='LA'
	       obslev(1:2)='QL'
	    else if (line(1:5).eq.'SIGMA' .and. line(16:17).eq.'52') then
	       obstyp(3:4)='PR'
	       obslev(3:4)='QL'
	    else if (line(1:5).eq.'SIGMA' .and. line(15:17).eq.'100') then
	       obstyp(5:6)='XO'
	       obslev(5:6)='QL'
	    else if (line(1:5).eq.'MBIAS' .and. line(15:17).eq.'900') then
	       obstyp(5:6)='XO'
	       obslev(5:6)='OP'
	    else if (line(1:5).eq.'MBIAS' .and. line(16:17).eq.'99') then
	       obstyp(3:4)='RA'
	       obslev(3:4)='OP'
	    else if (line(1:5).eq.'EARTH') then
	       read (line(25:44),*) gm
	    endif
         enddo
         if (buffer(1).ne.-8d9) goto 1310
      enddo

      if (qualit.gt.6) then
	 qualit=1
      else
	 qualit=0
      endif
      if (satid(isat).eq.9105001) then
	 satel='ERS1'
	 modid=5
      else if (satid(isat).eq.9502101) then
	 satel='ERS2'
	 modid=5
      else if (satid(isat).eq.0105501) then
         satel='JAS1'
         modid=5
      endif
      if (type.eq.4) then
	 prcrec='DSIDP '//satel//'.ORB.'//ctype//
     |		'   POSVEL DUT/DEOS Orbit'
	 write (6,682) prcrec
	 write (prcrec,681) 'STATE ',et0-et2000,et1-et2000,
     .		obstyp,obslev,modid,relid,0,0,0,qualit,
     .		mod(et0*secday,1800d0),' DUT/DEOS Orbit'
681   format (a6,2f7.1,2a6,2i2,3i4,i1,f6.3,a)
	    call convrec(prcrec)
	    write (6,682) prcrec
      endif
      xmjd0=0
      quali=2

* Process Data buffers

  100 read (30,err=1300,end=1330) buffer
      if (buffer(1).ne.9d9) then
      else if (type.eq.4) then
	 rewind (30)
	 type=5
	 goto 10
      else
	 goto 200
      endif

* Store date and start-time

      idate=int(buffer(2)/1d6)
      mjd=mdate(2,idate)
      time=buffer(2)-idate*1d6
      ihour=int(time/1d4)
      time=time-ihour*1d4
      imin=int(time/1d2)
      time=time-imin*1d2
      time=ihour*3600+imin*60+time
      if (jdyear.eq.0) then
	 jdyear=mdate(2,idate/10000*10000+0101)-30000
      endif

      ntb=nint(buffer(5))
      ibuf=2*ntimbf-1

      if (type.eq.2.and.mjd0.eq.0) then
	 mjd0=mjd
	 time0=time
	 call unlink("/tmp/out.pon")
	 open (10,file="/tmp/out.pon",status='new')
      endif
      if (xmjd0.eq.0d0) then
	 xmjd0=mjd+time/secday
	 xmjdp0=xmjd0+1d0-35d0/1002
      endif

      do idat=1,ntb
	 do i=1,6
	    xyz(i)=buffer(ibuf+ioff+i)
	 enddo
	 do i=1,3
	    geo(i)=buffer(ibuf+ioff+6+i)
	 enddo
	 do i=1,6
	    ecf(i)=buffer(ibuf+ioff+9+i)
	 enddo
         grhran=buffer(5+ntimbf+idat)/degrad
	 xpole=buffer(ibuf+ioff+16)
	 ypole=buffer(ibuf+ioff+17)
	 if (time.ge.secday) then
	    time=time-secday
	    mjd=mjd+1
	 endif
	 xmjd=mjd+time/secday
	 if (allarc) then
	    quali=0
	 else if (xmjd-xmjd0.gt.5.5d0) then
	    goto 100
	 else if (xmjd.gt.xmjdp0 .and. xyz(6).gt.0 .and. zdot.lt.0) then
	    quali=2-quali
	    xmjdp0=xmjd+3.5d0-35d0/1002
	 endif

	 et=(buffer(4)+buffer(5+idat))/secday
         ibuf=ibuf+nwdsat

* Determine equator flag

	 ascarc=0
	 if (xyz(3).gt.0 .and. z.lt.0) ascarc=1
	 z=xyz(3)
	 zdot=xyz(6)

         if (type.eq.1) then

* Write POE format

	    call convtime(mjd,time,iyear,imonth,iday,ihour,imin,sec)
	    write (6,640) iyear,imonth,iday,ihour,imin,sec,
     .		grhran,xpole,ypole,(et-jdyear)
	    write (6,650) xyz
	    write (6,650) ecf
            write (6,655) atc1,atc2

         else if (type.eq.2) then

* Write PON format

	    call convtime(mjd,time,iyear,imonth,iday,ihour,imin,sec)
            xpole=xpole*masrad
            ypole=ypole*masrad
            call rotate(1,-ypole,ecf(1),ecf(1))
            call rotate(2,-xpole,ecf(1),ecf(1))
            call rotate(1,-ypole,ecf(4),ecf(4))
            call rotate(2,-xpole,ecf(4),ecf(4))
            call xyzgeo(ecf(1),r,geo(1),geo(2),geo(3))
	    geo(1)=geo(1)/degrad
	    geo(2)=geo(2)/degrad
	    if (geo(2).lt.0) geo(2)=geo(2)+360
	    write (10,660) 1900+iyear,imonth,iday,ihour,imin,sec,
     .		geo(2),geo(1),geo(3),xyz,ecf

         else if (type.eq.3) then

* Write DUT format

            call rotate(1,-ypole*masrad,ecf(1),ecf(1))
            call rotate(2,-xpole*masrad,ecf(1),ecf(1))
            call xyzgeo(ecf(1),r,geo(1),geo(2),geo(3))
	    geo(1)=geo(1)/degrad
	    geo(2)=geo(2)/degrad
	    if (geo(2).lt.0) geo(2)=geo(2)+360
	    write (6,670) mjd,time,grhran,xpole,ypole,geo,xyz

* Write PRC format (STINER cards)

	 else if (type.eq.4) then
	    call vecele(xyz,elem,gm)
	    u=elem(5)+elem(6)
	    pitch=-0.335*sin(u)*cos(u)
	    roll = 0.050*sin(u)
	    yaw  = atan(0.0683*cos(u))/degrad
*    write (6,*) (xyz(i),i=1,3)
	    call j2000(-1,xmjd,0d0,0d0,6,xyz,ecf)
*    write (6,*) (ecf(i),i=1,3)
*    call j2001(xmjd,xyz,ecf,6,2)
*    write (6,*) (ecf(i),i=1,3)
	    iday=int(et)
	    if (quali.ne.2) then
	       write (prcrec,680) 'STINER',satid(isat),orbtyp,
     .		iday-et2000,(et-iday)*secday,ecf,
     .		roll,pitch,yaw,ascarc,0,quali,0
680   format (a6,i7,a1,f7.1,f12.6,3f13.3,3f12.6,3f7.3,i2,i3,i1,i4)
	       call convrec(prcrec)
	       write (6,682) prcrec
	    endif
682   format (a129)

* Write PRC format (STTERR cards)
* Change ECF from instantaneous pole to mean pole.

	 else if (type.eq.5) then
	    call vecele(xyz,elem,gm)
	    u=elem(5)+elem(6)
	    pitch=-0.335*sin(u)*cos(u)
	    roll = 0.050*sin(u)
	    yaw  = atan(0.0683*cos(u))/degrad
            call rotate(1,-(ypole-ypole0)*masrad,ecf(1),ecf(1))
            call rotate(2,-(xpole-xpole0)*masrad,ecf(1),ecf(1))
	    iday=int(et)
	    if (quali.ne.2) then
	       write (prcrec,680) 'STTERR',satid(isat),orbtyp,
     .		iday-et2000,(et-iday)*secday,ecf,
     .		roll,pitch,yaw,ascarc,0,quali,0
	       call convrec(prcrec)
	       write (6,682) prcrec
	    endif
         endif

         time=time+dtime
      enddo
      goto 100

200   continue
      if (type.eq.2) then
	 close (10)
	 call convtime(mjd0,time0,iyear,imonth,iday,ihour,imin,sec)
	 write (6,661) 1900+iyear,imonth,iday,ihour,imin,sec
	 call convtime(mjd,time-dtime,iyear,imonth,iday,ihour,imin,sec)
	 write (6,661) 1900+iyear,imonth,iday,ihour,imin,sec
	 write (6,550)
	 call system('cat /tmp/out.pon;/bin/rm -f /tmp/out.pon')
      endif
      goto 9999

* Formats

  550 format (a)
  640 format ('0.',i4.4,4i2.2,'0000D+12',5D22.16)
  650 format (6D22.16)
  655 format (22i1,4d22.16)
  660 format (1x,i4,4(1x,i2),1x,f10.7,2(1x,f8.3),13(1x,f14.4))
  661 format (1x,i4,4(1x,i2),1x,f10.7,$)
  670 format (i5,f9.3,f13.9,2f7.2,2f13.9,f13.4,3f14.4,3f12.6)

* Errors

*1300  call perror(program)
1300  write(0,*) program // ': error reading from ' // filenm
      goto 9999

1310  write (0,1311) program,'file is not a G2T file'
1311  format (a,': ',a)
      goto 9999

1320  write (0,1311) program,
     .'latitude, longitude or height not available'
      goto 9999

1330  write (0,1311) program,'premature end of file'
      goto 9999

1340  write (0,1311) program,'multiple satellites not allowed'
      goto 9999

1350  write (0,1311) program,'can''t open G2T file'
      goto 9999

1360  write (0,550) 'g2tpoe: unrecognized command: '//filenm
9999  end

      subroutine convtime(mjd,time,iy,mm,id,ih,im,sec)
      integer is,id,ih,im,iy,mm,mjd,mdate
      real*8 sec,time
      is=nint((mjd-46066)*86400d0+time)
      id=is/86400
      is=is-id*86400
      ih=is/3600
      is=is-ih*3600
      im=is/60
      is=is-im*60
      sec=is
      id=mdate(1,id+46066)
      iy=id/10000
      mm=id/100-iy*100
      id=mod(id,100)
      if(iy.gt.57) then 
        iy = iy + 1900
      else
        iy = iy + 2000
      endif
      end

      subroutine convrec(a)
      character*(*) a
      integer i,j,l,check

      check=0
      j=0
      do i=1,len(a)
	 if (a(i:i).ne.'.') then
	    j=j+1
	    a(j:j)=a(i:i)
	    l=ichar(a(j:j))
	    if (j.ge.21 .and. j.le.120 .and.
     .		l.ge.48 .and. l.le.57)  check=check+(l-48)
	 endif
      enddo
      if (a(1:6).eq.'STTERR' .or. a(1:6).eq.'STINER')
     .		write (a(121:123),'(i3)') check
      end
