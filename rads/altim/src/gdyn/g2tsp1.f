      program g2todr

* Convert GEODYN II Trajectory file (G2T) to SP1 format.
*
*  4-Jul-2002 - Created by Eelco Doornbos 
************************************************************************
      implicit none

      real*8 buffer(2048)
      character*80 filenm, arg
      integer*4 iarg, iargc, idat, i
      integer*4 na, nwdsat, ntimbf, dtime, ioff, ibuf

      real*8 sec85, ymdhms

      integer*4 ntb, mjdsbf

      real*8 pi, rad, rmul

* SP1 header records 

      integer*4 year, month, day, hour, minute
      real*8 second, interval, fracday, gpsmjd, gpsymd
      integer*4 mjd, nepoch, gpsweek 
      integer isat(34)
      data isat /34*0/

* SP1 data records

      real*8 ecf(6), xpole, ypole

      pi=4*datan(1.0d0)
      rad=pi/180.0d0
      rmul=rad/3600.0d3

* Scan the command line

      do iarg=1,iargc()
        call getarg(iarg,arg)
        filenm=arg
      enddo

* Open G2T file

      open (30,status='old',file=filenm,form='unformatted',err=1350)

* Read G2T header buffer

      rewind (30)
      read (30) buffer
      if (buffer(1).ne.-9d9) write (0,*) 'buffer(1)=',buffer(1)
      na=nint(buffer(2))
      if (buffer(7).ne.1) then
        stop 'More than one satellite in geodyn orbit file.'
      end if
      nwdsat=nint(buffer(8))
      ntimbf=nint(buffer(10))
      dtime=buffer(19)
      ioff=0
      do i=202,207
         if (buffer(i).gt.0) ioff=ioff+1
      enddo 
      if (buffer(210).le.0 .or. buffer(211).le.0 .or.
     |    buffer(212).le.0) then
        stop 'No ECF XYZ coordinates in geodyn orbit file.'
      end if
      if (buffer(213).le.0 .or. buffer(214).le.0 .or.
     |    buffer(215).le.0) then
        stop 'No ECF XYZ velocities in geodyn orbit file.'
      end if
      isat(1)=nint(buffer(301))

* Write SP1 header

      gpsmjd = (buffer(15) + buffer(16) - 32.184d0 - 19.0d0) / 
     | 86400.0D0 + 30.0D3
      gpsymd = ymdhms(sec85(1,gpsmjd))
      mjd = int(gpsmjd)
      gpsweek = int( (dble(mjd) - 44244.0D0) / 7.0D0 )

      call ymdsplit(gpsymd,year,month,day,hour,minute,second)
      fracday = (hour*3600.0D0+minute*60.0D0 + second) / 86400.0d0

      interval = buffer(19)
      nepoch = buffer(20)

      write (*,910) ' # ', year, month, day, hour, minute, second, 
     | interval, mjd, fracday, nepoch, 'F', 'DEOS'
      write (*,920) ' + ', 1, isat, 0, gpsweek

 910  format (a3,1x,i4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,f10.7,1x,
     | f14.7,1x,i5,1x,f15.13,1x,i6,a1,a4)
 920  format (a3,i2,1x,34i2,i2,i4)

* Read Alphanumeric header(s)

      read (30,err=1330,end=1330) buffer
      if (buffer(1).ne.-8d9) stop 'File is not a geodyn orbit file'
      do i=2,na
         read (30,err=1330,end=1330) buffer
      enddo

  10  read (30,err=1330,end=1330) buffer
      if (buffer(1).eq.9d9) goto 9999
        mjdsbf = buffer(4)
        ntb = nint(buffer(5))
        ibuf = 2 * ntimbf - 1

        do idat=1, ntb
          gpsmjd = (mjdsbf + buffer(5+idat) - 32.184d0 - 19.0d0) /
     |              86400.0D0 + 30.0D3
          gpsymd = ymdhms(sec85(1,gpsmjd))
          call ymdsplit(gpsymd,year,month,day,hour,minute,second)

          do i=1,6
            ecf(i) = buffer(ibuf+ioff+9+i) / 1d3
          enddo 
          xpole = buffer(ibuf+ioff+16) * rmul
          ypole = buffer(ibuf+ioff+17) * rmul
          call rotate(2,-xpole,ecf,ecf)
          call rotate(1,-ypole,ecf,ecf)

          ibuf = ibuf + nwdsat

          write (*,930) ' * ', year, month, day, hour, minute, second
          write (*,940) ' SV', isat(1), ecf

 930      format (a3,1x,i4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,f10.7)
 940      format (a3,i2,3f13.6,3f12.8)

        enddo

      goto 10

 1330 continue
      stop 'Unable to read from GEODYN orbit file'

 1350 continue 
      stop 'Unable to open GEODYN orbit file'

 9999 continue

      write (*,'(a3)') 'EOF'

      end

********************************************************************************

      subroutine ymdsplit(ymd,year,month,day,hour,minute,second)
      real*8 ymd, second
      integer*4 year,yy,month,day,hour,minute

      yy = int(ymd/1e10)
      if (yy.le.56) then
        year = yy + 2000
      elseif (yy.gt.57.and.yy.le.1957) then
        year = yy + 1900
      else
        year = yy
      endif
      month = int((ymd-yy*1d10)/1d8)
      day = int((ymd-yy*1d10-month*1d8)/1d6)

      hour = int((ymd-yy*1d10-month*1d8-day*1d6)/1d4)
      minute = int((ymd-yy*1d10-month*1d8-day*1d6-hour*1d4)
     |         / 1d2)
      second = ymd-yy*1d10-month*1d8-day*1d6-hour*1d4 -
     | minute*1d2
 
      end
