C---*-$--1----*----2----*----3----*----4----*----5----*----6----*----7-<--*----8
      program dblepass

C This program deletes passes from the file gbfdata2, which occure both in the
C files gbfdata1 and gbfdata2. The program starts with making an inventory off
C all relevant passes in the file gbfdata1. For the begin and end time of the
C passes the times of the first and last measurement are taken, substracted and
C added by the value deltat. Then each measurement of gbfdata2 is read and
C checked if it coincides with one of the passes of gbfdata1, if not this
C measurement is written to the output file gbfdata2.new.
C
C input files:
C
C gbfdata1  (if10)  - File with measurements in gbf format.
C
C gbfdata2  (if11)  - File with measurements in gbf format, from which data is
C                     deleted which coincide with data passes in gbfdata1.
C
C output files:
C
C standard-output (terminal) (if6)
C
C gbfdata2.new (if10) - File with remaining measurements from gbfdata2.
C
C This program has an optional argument:
C
C      deltat=value - where value is the number of minutes the passes are
C                     extended with. Default value is 50.0.
C
C Created by G.J.Mets, may 1995.
C

      implicit double precision ( a - h , o - z )

      parameter ( maxpas = 20000)  

      integer*2 mtype 
      character *30 carg
      logical lexist
      dimension idata(17)

      equivalence (idata(1),idsat),(idata(2),mtype),(idata(3),istat),
     .            (idata(5),mjd),(idata(6),tfrac)

      common / passes / t1pas(maxpas), t2pas(maxpas), mtppas(maxpas), 
     .                  isapas(maxpas), istpas(maxpas)
      common / files / if6, if10, if11, if20

      data if6 , if10, if11, if20 / 6, 10, 11, 20 /
      data deltat / 50.0d0 /

c Open input and output files.

      inquire ( file = 'gbfdata1' , exist = lexist )
      if ( lexist ) then
        open ( if10, file = 'gbfdata1', form='unformatted' )
      else
        write ( if6, * ) 'ERROR: File gbfdata1 does not exist.'
        write ( if6, * ) '       Execution stops in program dblepass.'
        stop 99
      endif
      inquire ( file = 'gbfdata2' , exist = lexist )
      if ( lexist ) then
        open ( if11, file = 'gbfdata2', form='unformatted' )
      else
        write ( if6, * ) 'ERROR: File gbfdata2 does not exist.'
        write ( if6, * ) '       Execution stops in program dblepass.'
        stop 99
      endif
      open ( if20, file = 'gbfdata2.new' , form='unformatted' )

c Read argument if available.

      narg = 0
      narg = iargc()
      if ( narg .ne. 0 ) then
        call getarg ( 1, carg )
        if ( carg(1:7) .ne. 'detat=' ) then
          write ( if6, * ) 'ERROR: Unknown argument: ', carg
          write ( if6, * ) '       Execution stops in program dblepass.'
          stop 99
        endif
        read ( carg(8:30), * ) deltat
      endif

c Convert deltat from minutes to fraction of day.

      deltat = deltat / 1440.0d0

c Make an inventory of the passes in the file gbfdata1.

      call invent ( deltat, nrpass )

c Start reading the data from gbfdata2.

      nrdel = 0
      istart = 1

 10   read ( if11, end = 100 ) idata

      tmeas = dble ( mjd ) + tfrac

c Check if measurement coincides with one of the passes from gbfdata1.
c If not write it to gbfdata2.new.

      do i = istart, nrpass
        if ( tmeas .gt. t2pas(i) ) then
          istart = istart + 1
        else
          if ( tmeas .gt. t1pas(i) ) then
            if ( istat .eq. istpas(i) ) then
              if ( mtype .eq. mtppas(i) ) then
                if ( idsat .eq. isapas(i) ) then
                  nrdel = nrdel + 1
                  goto 10
                endif
              endif
            endif
          else
            if ( tmeas + 3.0d0 * deltat .lt. t1pas(i) ) goto 50
          endif
        endif
      enddo

 50   write ( if20 ) idata

      goto 10

 100  close ( if10 )
      close ( if11 )
      close ( if20 )
      write ( if6, * ) nrdel, ' measurements were deleted.'
      write ( if6, * ) 'Normal end of program dblepass.'

      end

C---*-$--1----*----2----*----3----*----4----*----5----*----6----*----7-<--*----8

      subroutine invent ( deltat, nrpass )

C This subroutine makes an inventory of all passes in the file gbfdata1. The
C inventory is started after the time of the first measurement in the file
C gbfdata2 and is ended after the time of the last measurement. After the
C inventory the passes are sorted on time taking the end time (t2pass) of the
C passes.

      implicit double precision ( a - h , o - z )

      parameter ( maxpas = 20000)

      integer * 2 mtype , indtim

      common / passes / t1pas(maxpas), t2pas(maxpas), mtppas(maxpas), 
     .                  isapas(maxpas), istpas(maxpas)
      common / files / if6, if10, if11, if20

c Read first observation from gbfdata2.

      read ( if11 ) idsat, mtype, indtim, istat, iprepr, mjd, fract

c Determine start time of the inventory.

      timbeg = dble ( mjd ) + fract - deltat

c Read last observation from gbfdata2.

 5    read ( if11, end = 6 ) idsat, mtype, indtim, istat, iprepr,
     .                       mjd, fract
      goto 5

c Determine end time of the inventory.

 6    timend = dble ( mjd ) + fract + deltat
      rewind if11

c Start the inventory.

      nrpass = 0

  10  read ( if10, end = 100 ) idsat, mtype, indtim, istat,
     .                         iprepr, mjd, fract
      time = dble ( mjd ) + fract

      if ( time .lt. timbeg ) goto 10
      if ( time .gt. timend ) goto 100

c  check whether taken during "existing" pass

      do 20 i = nrpass , 1 , -1
        if ( istat .eq. istpas(i) ) then
          if ( idsat .eq. isapas(i) ) then
            if ( mtype .eq. mtppas(i) ) then
              dt = time - t2pas(i)
              if ( dt .le. deltat ) then
                t2pas(i) = time
                goto 10
              end if
            end if
          end if
        end if
  20    continue

c  new pass

      call plus1 ( if6, nrpass , maxpas , 'ERROR' , 'maxpas' )
      t1pas(nrpass)  = time
      t2pas(nrpass)  = time
      mtppas(nrpass) = mtype
      isapas(nrpass) = idsat
      istpas(nrpass) = istat
      goto 10

 100  continue

c extend all passes with deltat on both ends.

      do i = 1, nrpass
        t1pas(i) = t1pas(i) - deltat
        t2pas(i) = t2pas(i) + deltat
      enddo

c Sort the passes on time, taking t2pas

      do i = 2, nrpass
        ind = i
        if ( t2pas(ind) .lt. t2pas(ind-1) ) then
          t2hlp  = t2pas(ind)
          t1hlp  = t1pas(ind)
          mtphlp = mtppas(ind)
          isahlp = isapas(ind)
          isthlp = istpas(ind)
 110      t2pas(ind)  = t2pas(ind-1)
          t1pas(ind)  = t1pas(ind-1)
          mtppas(ind) = mtppas(ind-1)
          isapas(ind) = isapas(ind-1)
          istpas(ind) = istpas(ind-1)
          ind = ind - 1
          if ( ind .eq. 1 ) goto 120
          if ( t2hlp .lt. t2pas(ind-1) ) goto 110
 120      t2pas(ind)  = t2hlp
          t1pas(ind)  = t1hlp
          mtppas(ind) = mtphlp
          isapas(ind) = isahlp
          istpas(ind) = isthlp
        endif
      enddo

      end
