      program dorisedit
  
      implicit none

      character*5  currstat

      integer*4    ipass, imeas, maxmeas
      parameter    (maxmeas=250)
      character*96 line, l(maxmeas)
      
      integer*4    npre, npost, nfedit, nftropo, nfiono, ntime
      integer*4    nout, ninput, ngap, nshort
      common /lerror/ npre, npost, nfedit, nftropo, nfiono, ntime, 
     |                nout, ninput, ngap, nshort

      real*8 timetag, ptimetag

* Initialize

      imeas = 0
      ipass = 0
      ninput = 0
      nout = 0
      npre = 0
      npost = 0
      nfedit = 0
      nftropo = 0
      nfiono = 0
      ntime = 0
      nshort = 0
      ngap = 0
      currstat = 'xxxxx'

* Read data lines
 
  10  read (*,'(a96)',end=99) line

* Skip header lines (for Envisat DOR_DOP files)

      if (line(8:9).ne.'39') goto 10

      ninput = ninput + 1

* Start a new pass if station changes or if there is a gap of more than 20 minutes

      if (line(12:16).eq.currstat .and.
     |    timetag(line).lt.ptimetag+20.0D0*60.0D0)  then
        imeas = imeas + 1
        l(imeas) = line
      else
        if (currstat.ne.'xxxxx') call dumppass(l,imeas)
        ipass = ipass + 1
        imeas = 1
        l(imeas) = line
      endif

      currstat = line(12:16)
      ptimetag = timetag(line)

      goto 10

  99  continue

      if (currstat.ne.'xxxxx') call dumppass(l,imeas)

      write(7,200) 'lines read:', ninput
      write(7,200) 'points edited during pre-processing:', npre
      write(7,200) 'points edited during post-processing:', npost
      write(7,200) 'points edited due to time-tag overlaps:', ntime
      write(7,200) 'points edited due to short pass length:', nshort
      write(7,200) 'points edited due to large gap in pass:', ngap
      write(7,201) 'lines written:', nout, 
     |              '(',int((dble(nout)/dble(ninput))*100.0D0),'%)'
      write(7,200) 'number of passes:', ipass
 200  format (a40,1x,i6)
 201  format (a40,1x,i6,1x,a1,i3,a2)

      end


********************************************************************************

      logical function checkline(line)

* DORIS format elements

      character*96 line
      character*7  satnm
      integer*4    tsys1, tsys2
      character*5  stationid
      integer*4    prepro1
      integer*4    prepro2
      integer*4    prepro3
      integer*4    count

      integer*4    npre, npost, nfedit, nftropo, nfiono, ntime
      integer*4    nout, ninput, ngap, nshort
      common /lerror/ npre, npost, nfedit, nftropo, nfiono, ntime, 
     |                nout, ninput, ngap, nshort

      checkline = .true.

      satnm = line(1:7)
      read(line(10:10),*,err=19) tsys1
      read(line(11:11),*,err=19) tsys2
      stationid = line(12:16) 
      read(line(33:33),*,err=19) prepro1
      read(line(34:34),*,err=19) prepro2
      read(line(35:35),*,err=19) prepro3
      read(line(36:45),*,err=19) count

* Check if iono correction applied
      if (prepro1.eq.1) then
      elseif (prepro1.gt.1.or.prepro1.lt.0) then
        checkline = .false.
      endif


* Check if tropo correction applied
      if (prepro2.eq.1) then
      elseif (prepro2.gt.1.or.prepro2.lt.0) then
        checkline = .false.
      endif

* Check pre- and post-processing editing flag
      if (prepro3.eq.1) then
        checkline = .false.
        npre = npre + 1
      elseif (prepro2.eq.2) then
        checkline = .false.
        npost = npost + 1
      elseif (prepro2.gt.2.or.prepro2.lt.0) then
        checkline = .false.
        nfedit = nfedit + 1
      endif

      return

 19   print *, 'Error reading DORIS data from the following line:'
      print *, line
      stop

      end

********************************************************************************

      real*8 function timetag(line)
      character*96 line
      integer*4 iyear, idoy, isecofday, secfrac, count
      real*8 yydoyfr, sec85

      read(line(17:18),*,err=19) iyear
      read(line(19:21),*,err=19) idoy
      read(line(22:26),*,err=19) isecofday
      read(line(27:32),*,err=19) secfrac
      read(line(36:45),*,err=19) count      
      yydoyfr = iyear * 1.0d3 + idoy + isecofday / 86400d0 +
     |          secfrac / 1.0d6
      timetag = sec85(3,yydoyfr) + count / 1.0d7

      return

 19   print *, 'Error reading DORIS data from the following line:'
      print *, line
      stop

      end 

********************************************************************************

      subroutine dumppass(lines,n)

      integer*4    n, i, j
      character*96 lines(*)
      real*8 timetag
      logical checkline

      integer*4    npre, npost, nfedit, nftropo, nfiono, ntime
      integer*4    nout, ninput, ngap, nshort
      common /lerror/ npre, npost, nfedit, nftropo, nfiono, ntime, 
     |                nout, ninput, ngap, nshort

* Check if pass longer than 4 minutes

      if (timetag(lines(n)) - timetag(lines(1)) .lt. 240.0d0) then
        nshort = nshort + n
        return
      endif

* Check if pass does not contain gaps of over 4 minutes

      do i=2,n
        if(timetag(lines(i))-timetag(lines(i-1)).gt.240.0D0) then
          ngap = ngap + n
          return
        endif
      enddo

* Start output (always start with first line of pass)

      nout = nout + 1
      j = 1
      if(checkline(lines(1))) write(*,'(a96)') lines(1)
 
      do i=2, n
        if(timetag(lines(i))-timetag(lines(j)).lt.6.0D0) then
          ntime = ntime + 1
        else
          if(checkline(lines(i))) then
            nout = nout + 1
            j = i
            write(*,'(a96)') lines(i)
          endif
        endif
      enddo

      return

      end
