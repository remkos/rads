      program mrt2gbf
c*********************************************************************
c
c     This program converts laser ranging measurements from the
c     MERIT II format to the geodyn binary format. in addition to
c     the strict conversion, some new information may be added:
c     * if the time scale of the input data records is in utc(bih),
c       (code "7"), the time-tag is corrected to the utc(usno) scale
c       (code "3") by substracting an interpolated value for utc(bih)
c       - utc(usno). values for the latter are taken from a tso-file
c       with the appropriate information.
c     * the measurements are corrected for the center-of-mass offset,
c       if this correction has not been applied to the input records
c       yet (lageos, starlette, ajisaj, ers-1, and t/p only).
c     * the field for the tropospheric refraction correction is filled
c       with this value if it is available from the input data record.
c       otherwise, the appropriate value for the laser wavelength is
c       encoded in this field.
c
c     input files:
c     stdin- a dataset with (laser ranging) measurements in the
c            MERIT II format.
c     ft11 - a dataset with information on the utc(bih) - utc(usno)
c            time difference.
c     ft12 - a dataset with information on the wavelength of the laser
c            systems.
c
c     output file:
c     ft20 - a dataset with translated laser range measurements in the
c            geodyn binary format.
c
c     r. noomen, june 23, 1989.
c     r. noomen, may 2, 1990.     (ajisaj added)
c     r. noomen, july 25, 1991.   (ers-1 added; cmass correction = 0 ! )
c     r. noomen, september 23, 1991.   (cmass correction ers-1 = 0.061 m
c                                       zero value in data field; cmass
c                                       indicator "on" for offset)
c     r. noomen, june 12, 1992.  (adjusted for 251 mm cmass correction
c                                 for lageos1)
c     r. noomen, october 30, 1992. (lageos-2 added)
c     r. scharroo, 18 nov 1993, ported to convex.
c
c*********************************************************************
c
      implicit double precision ( a - h , o - z )
      integer   * 2 j5
      real          std,trop,iono,displ,cmas
      character * 1 c1
      character * 80 filenm

      parameter ( nsat = 6 )
      parameter ( ndim = 1000 )

      dimension idata(17) , j5(2)
      dimension isat(nsat) , cmass(nsat) , imass(nsat)
      dimension itwo(31) , dmjdy(99)

      equivalence
     .  ( idata(1) , j1 ) ,    ( idata(2) , j5(1) ) ,
     .  ( idata(3) , j9 ) ,    ( idata(4) , j13 ) ,
     .  ( idata(5) , j17 ) ,   ( idata(6) , tfrac ) ,
     .  ( idata(8) , range ) , ( idata(10) , j37 ) ,
     .  ( idata(11) , j41 ) ,  ( idata(12) , std ) ,
     .  ( idata(13) , trop ) , ( idata(14) , j53 ) ,
     .  ( idata(15) , iono ) , ( idata(16) , displ ) ,
     .  ( idata(17) , cmas )

      common / meas / ntot(8) , nmeas(8,9999)
      common / utc  / nutc , xutc(ndim) , tutc(ndim)
      common / wav  / nwav , istwav(ndim) ,
     .                t1wav(ndim) , t2wav(ndim) , iwav(ndim)

      data ifin , ifout / 5 , 20 /
      data icur , vlight / 1 , 299792458.0d0 /
      data isat , cmass , imass /
     .     7603901 , 9207002 , 7501001 , 8606101 , 9105001 , 9205201 ,
     .     0.251d0 , 0.251d0 , 0.075d0 , 1.01d0  , 0.061d0 , 0.000d0 ,
     .     0       , 0       , 0       , 0       , 1       , 1 /

      call getarg (1,filenm)
      if (filenm.eq.' ') then
	 write (6,600)
600      format ('usage: mrt2gbf gbf'//
     .     'converts MERIT II data set on standard input to GBF dataset')
	 stop
      endif
      open (ifout,file=filenm,status='new',form='unformatted')

      call initdt ( 'times.data' )
      call initwl ( 'system.data' )
      call inityr ( dmjdy )
      call init2 ( 31 , itwo )

      j5(1) = 20
      j37   = 0
      j41   = 0
      iono  = 0.0
      displ = 0.0

  10  read   ( ifin , 20 , end = 100 ) i1 , i8 , i10 , i13a , i13b ,
     .         i25 , i29 , i31 , i33 , i40 , i46a , i46b ,
     .         i58 , i65 , i69 , i74 , i78 , i81 , i86 ,
     .         i92 , i97 , i105 , i111 , i115 , i116 ,
     .         i120 , i121 , i122 , i123 , i124 , i125 , i126 ,
     .         i127 , i128 , i129 , c1
  20  format ( 1i7 , 1i2 , 1i3 , 2i6 , 1i4 , 1i2 , 1i2 , 1i7 , 1i6 ,
     .         2i6 , 1i7 , 1i4 , 1i5 , 1i4 , 1i3 , 1i5 , 1i6 , 1i5 ,
     .         1i8 , 1i6 , 1i4 , 1i1 , 1i4 , 10i1 , 1a1 )

      if ( i120 .eq. 3 ) then
        nmeas(2,i25) = nmeas(2,i25) + 1
        ntot(2)      = ntot(2)      + 1
        goto 10
      end if
      if ( i121 .gt. 7 ) then
        nmeas(3,i25) = nmeas(3,i25) + 1
        ntot(3)      = ntot(3)      + 1
        goto 10
      end if
      if ( i8 .lt. 57 ) then
        nmeas(4,i25) = nmeas(4,i25) + 1
        ntot(4)      = ntot(4)      + 1
        goto 10
      end if

      date = dmjdy(i8) + dble ( i10 ) + dble ( i13a ) / 86400.0d1
     .       + dble ( i13b ) / 86400.0d7
      if ( i121 .eq. 7 ) then
        deltat = difft ( date , icur , i25 )
        date   = date - deltat / 86400.0d6
        i121=3
      end if

      j1    = i1
      j5(2) = 10 * i120 + i121
      j9    = i25
      if ( i69 .ne. 0 ) then
        if ( i123 .eq. 0 ) then
          itest = 4
        else
          itest = 5
        end if
      else
        if ( i123 .eq. 0 ) then
          itest = 0
        else
          itest = 1
        end if
      end if

      index = 0
      do 30 i = 1 , nsat
        if ( i1 .eq. isat(i) ) then
          index = i
          goto 31
        end if
  30  continue
      write ( 6 , * ) ' unknown satellite id. stop.'
      stop 99
  31  continue

      j13   = imass(index) * itwo(28) + 3 * itwo(25) + itest * itwo(19)
      j17   = date
      tfrac = date - dble ( j17 )
      range = ( dble ( i46a ) * 1.0d6 + dble ( i46b ) )
     .        * 0.5d-12 * vlight

c  new since june 12, 1992: test for correct cmass for lageos
      if ( index .eq. 1 .and. i124 .eq. 0 ) then
        if ( 1000 .lt. i86 .and. i86 .lt. 2000 ) then
          range = range - dble ( i86 ) * 1.5d-4 + cmass(index)
        else
          write ( 6 , * ) 'main: suspected cmass value for ' ,
     .                    'station ' , i25 , ' (continue)'
        end if
        cmas = cmass(index)
      end if

      if ( i124 .eq. 1 ) then
        range = range + cmass(index)
        if ( isat(index) .ne. 9105001 ) then
          cmas = cmass(index)
        else
          cmas = 0.0d0
        end if
c       nmeas(7,i25) = nmeas(7,i25) + 1
c       ntot(7)      = ntot(7)      + 1
      end if

      std = dble ( i58 ) * 0.5d-12 * vlight
      if ( i81 .ne. 0 ) then
        trop = dble ( i81 ) * 0.5d-12 * vlight
      else
        if ( i65 .eq. 0 ) then
          call wave(i65,i25,date,iwvcur,istcur,t1cur,t2cur)
          trop = ( dble ( i65 ) * 1.0d-3 + 99.0d0 ) * 1.0d-30
        else
          trop = ( dble ( i65 ) * 1.0d-4 + 99.0d0 ) * 1.0d-30
        end if
      end if
      ih  = i78
      it  = ( i74 + 5 ) / 10
      ip  = ( i69 + 5 ) / 10
      j53 = ih * itwo(24) + it * itwo(12) + ip
      write (ifout) idata
      nmeas(1,i25) = nmeas(1,i25) + 1
      ntot(1)      = ntot(1)      + 1
      goto 10

 100  call print
      stop

      end
c*********************************************************************
c
      double precision function difft ( date , icur , i25 )
c
      implicit double precision ( a - h , o - z )
c
      parameter ( ndim = 1000 )
c
      common / meas / ntot(8) , nmeas(8,9999)
      common / utc  / nutc , xutc(ndim) , tutc(ndim)
c
      if ( date .ge. tutc(icur) .and.
     .     date .le. tutc(icur+1) ) goto 100
      icur = icur + 1
      if ( icur .ne. nutc ) then
        if ( date .ge. tutc(icur) .and.
     .       date .le. tutc(icur+1) ) goto 100
      end if
      if ( date .lt. tutc(1) ) then
        nmeas(5,i25) = nmeas(5,i25) + 1
        ntot(5)      = ntot(5) + 1
        difft        = xutc(1)
        icur         = 1
        return
      end if
      if ( date .gt. tutc(nutc) ) then
        nmeas(6,i25) = nmeas(6,i25) + 1
        ntot(6)      = ntot(6) + 1
        difft        = xutc(nutc)
        icur         = nutc - 1
        return
      end if
c
      do 10 i = 1 , nutc - 1
        if ( date .ge. tutc(i) .and. date .le. tutc(i+1) ) then
          icur = i
          goto 100
        end if
  10  continue
c
 100  difft = xutc(icur) + ( xutc(icur+1) - xutc(icur) )
     .        * ( date - tutc(icur) ) / ( tutc(icur+1) - tutc(icur) )
      if ( icur + 1 .gt. nutc ) icur = icur - 1
      return
c
      end
c*********************************************************************
c
      subroutine init2 ( n , itwo )
c
      implicit double precision ( a - h , o - z )
      dimension itwo(n)
c
      itwo(1) = 2
      do 10 i = 2 , n
        itwo(i) = 2 * itwo(i-1)
  10  continue
      return
c
      end
c*********************************************************************
c
      subroutine initdt ( filenm )
c
      implicit double precision ( a - h , o - z )
c
      parameter ( ndim = 1000 , ifile = 10 )
      logical exist
      dimension i(2) , x(2)
      character*(*) filenm
c
      common / utc  / nutc , xutc(ndim) , tutc(ndim)
c
      inquire (file=filenm,exist=exist)
      if (exist) then
         open (ifile,file=filenm,form='formatted',status='old')
      else
         open (ifile,file='/user3/vlruweg/data/'//filenm,
     .   form='formatted',status='old')
      endif
      rewind (ifile)
      nutc = 0
      do 20 j = 1 , 3
        read   ( ifile , 10 , end = 999 )
  10    format ( a )
  20  continue
c
  30  read   ( ifile , * , end = 100 ) i , x
      if ( dabs ( x(1) ) .lt. 1.0d-10 .and.
     .     dabs ( x(2) ) .lt. 1.0d-10 ) goto 30
      call plus1 ( 6 , nutc , ndim , 'initdt' , 'tutc' )
      tutc(nutc) = dble ( i(2) )
      xutc(nutc) = x(1)
      goto 30
 100  continue
      close (ifile)
      return
c
 999  write ( 6 , * ) 'initdt: incorrect number of header records ',
     .  'in file with with utc(bih)-utc(usno) data. stop.'
      stop 99
c
      end
c*********************************************************************
c
      subroutine initwl ( filenm )
c
      implicit double precision ( a - h , o - z )
c
      parameter ( ndim = 1000 , ifile=10 )
      logical exist
      character*(*) filenm
c
      common / wav  / nwav , istwav(ndim) ,
     .                t1wav(ndim) , t2wav(ndim) , iwav(ndim)
c
      inquire (file=filenm,exist=exist)
      if (exist) then
         open (ifile,file=filenm,form='formatted',status='old')
      else
         open (ifile,file='/user3/vlruweg/data/'//filenm,
     .   form='formatted',status='old')
      endif
      rewind (ifile)
c
      i = 1
  10  read   ( ifile , 20 , end = 100 ) ihulp1 , ihulp2 ,
     .                                  iymd1 , iymd2 , iwav(i)
  20  format ( 1i4 , 4x , 1i1 , 31x , 2 ( 2x , 1i6 ) , 1i5 )
      istwav(i) = 10000 * ihulp2 + ihulp1
      call mjdate ( 2 , mjd1 , iymd1 , iflut1 , iflut2 , iflut3 )
      call mjdate ( 2 , mjd2 , iymd2 , iflut1 , iflut2 , iflut3 )
      t1wav(i) = dble ( mjd1 )
      t2wav(i) = dble ( mjd2 )
      call plus1 ( 6 , i , ndim , 'initwl' , 't1wav' )
      goto 10
c
 100  nwav = i - 1
      close (ifile)
      return
c
      end
c*********************************************************************
c
      subroutine inityr ( dmjdy )
c
      implicit double precision ( a - h , o - z )
      dimension dmjdy(99)
c
      do 10 i = 57 , 99
        call mjdate ( 2 , mjd , 10000 * i + 101 , if1 , if2 , if3 )
        dmjdy(i) = dble ( mjd ) - 1.0d0
  10  continue
      return
c
      end
c*********************************************************************
c
      subroutine print
c
      implicit double precision ( a - h , o - z )
c
      common / meas / ntot(8) , nmeas(8,9999)
c
      write  ( 6 , 110 )
 110  format ( ' overview of translated measurements:' , // ,
     .         ' station     measurements' , / , 1x , 24 ( '-' ) , / )
      do 130 i = 1 , 9999
        if ( nmeas(1,i) .ne. 0 ) write  ( 6 , 120 ) i , nmeas(1,i)
 120    format ( 1i6 , 5x , 1i10 )
 130  continue
      write  ( 6 , 140 ) ntot(1)
 140  format ( 11x , 14 ( '-' ) , / , 11x , 1i10 )
c
      if ( ntot(5) .ne. 0 ) then
        write  ( 6 , 520 ) 'before the first'
 520    format (///' warning: overview of measurements corrected' ,
     .          ' for the difference between utc(bih) and utc(usno),'/
     .          ' dated ',a,' value available from' ,
     .          ' the data file (this value was used):' //
     .          ' station     measurements' / 1x , 24('-') /)
        do 530 i = 1 , 9999
          if ( nmeas(5,i) .ne. 0 ) write  ( 6 , 120 ) i , nmeas(5,i)
 530    continue
      end if
c
      if ( ntot(6) .ne. 0 ) then
        write  ( 6 , 520 ) 'after the final'
        do 630 i = 1 , 9999
          if ( nmeas(6,i) .ne. 0 ) write  ( 6 , 120 ) i , nmeas(6,i)
 630    continue
      end if
c
      if ( ntot(7) .ne. 0 ) then
        write  ( 6 , 720 )
 720    format ( /// , ' warning: overview of measurements not yet' ,
     .           ' corrected for center-of-mass offset' ,
     .           ' (unknown satellite id.)' , / ,
     .           ' (geodyn can not apply this correction!):' , // ,
     .           ' station     measurements' ,
     .           / , 1x , 24 ( '-' ) , / )
        do 730 i = 1 , 9999
          if ( nmeas(7,i) .ne. 0 ) write  ( 6 , 120 ) i , nmeas(7,i)
 730    continue
      end if
c
      if ( ntot(8) .ne. 0 ) then
        write  ( 6 , 820 )
 820    format ( /// , ' warning: no proper wavelength could be' ,
     .           ' found for the following measurements' , / ,
     .           ' (a value of 600 nm was encoded in the translated' ,
     .           / , ' measurements):' , // ,
     .           ' station     measurements' ,
     .           / , 1x , 24 ( '-' ) , / )
        do 830 i = 1 , 9999
          if ( nmeas(8,i) .ne. 0 ) write  ( 6 , 120 ) i , nmeas(8,i)
 830    continue
      end if
c
      if ( ntot(2) .ne. 0 ) then
        write  ( 6 , 220 )
 220    format ( /// , ' overview of measurements deleted because of' ,
     .           ' time event references undefined in geodyn format:' ,
     .           // , ' station     measurements' ,
     .           / , 1x , 24 ( '-' ) , / )
        do 230 i = 1 , 9999
          if ( nmeas(2,i) .ne. 0 ) write  ( 6 , 120 ) i , nmeas(2,i)
 230    continue
      end if
c
      if ( ntot(3) .ne. 0 ) then
        write  ( 6 , 320 )
 320    format ( /// , ' overview of measurements deleted because of' ,
     .           ' epoch time scale undefined in geodyn format:' , // ,
     .           ' station     measurements' ,
     .           / , 1x , 24 ( '-' ) , / )
        do 330 i = 1 , 9999
          if ( nmeas(3,i) .ne. 0 ) write  ( 6 , 120 ) i , nmeas(3,i)
 330    continue
      end if
c
      if ( ntot(4) .ne. 0 ) then
        write  ( 6 , 420 )
 420    format ( /// , ' overview of measurements deleted because ' ,
     .           'of incorrect year of century:' , // ,
     .           ' station     measurements' ,
     .           / , 1x , 24 ( '-' ) , / )
        do 430 i = 1 , 9999
          if ( nmeas(4,i) .ne. 0 ) write  ( 6 , 120 ) i , nmeas(4,i)
 430    continue
      end if
      return
c
      end
c*********************************************************************
c
      subroutine wave ( inew , istnew , tnew , iold , istold ,
     .                  t1old , t2old )
c
      implicit double precision ( a - h , o - z )
c
      parameter ( ndim = 1000 )
c
      common / meas / ntot(8) , nmeas(8,9999)
      common / wav  / nwav , istwav(ndim) ,
     .                t1wav(ndim) , t2wav(ndim) , iwav(ndim)
c
      if ( istnew .eq. istold .and. tnew .ge. t1old
     .  .and. tnew .le. t2old ) then
        inew = iold
        return
      end if
c
      do 10 i = 1 , nwav
        if ( istwav(i) .ne. istnew ) goto 10
        if ( tnew .le. t1wav(i) .or. tnew .ge. t2wav(i) ) goto 10
        inew   = iwav(i)
        iold   = iwav(i)
        istold = istwav(i)
        t1old  = t1wav(i)
        t2old  = t2wav(i)
        return
  10  continue
c
      inew   = 600
      iold   = 0
      istold = 0
      t1old  = 0.0d0
      t2old  = 0.0d0
      nmeas(8,istnew) = nmeas(8,istnew) + 1
      ntot(8) = ntot(8) + 1
      return
c
      end
