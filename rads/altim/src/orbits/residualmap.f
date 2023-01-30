      program trackingmap

* Plots the availability of tracking data as a function of time and on
* a map of the world.
* Reads ODR files and (part of) lres output files.
* -
* 17 Jul 2002 - ED: First version

      implicit none

      integer maxpas, maxsta, maxplot, nplot, maxpoints
      parameter (maxpas=80000,maxsta=400,maxplot=100000,maxpoints=10000)
      character*80 staposfile, arg, device
      character*80 lresfile, odrfile, sat, cyc, title
      character*1 type
      integer i,j,k,iargc
      integer*4 jtrack, jtotal
      integer mjdref, statid(maxsta), station(maxpas), npas, nstapos
      integer splotted(maxsta),nplotted
      integer*4 nsta, nslr, ndor, npra, npdop, npran
      real*8  pos(6,maxsta), sig(6,maxsta)
      real*8  lat, lon, height, r, xyz(3), alt
      real*8  sec1(maxpas), sec2(maxpas), secb1, secb2, dt, mjd
      real*8  mean(maxpas), rms(maxpas)
      real*8  splot(maxplot)
      integer*4 nstaplot(maxplot), npoints
      real x(maxpoints), y(maxpoints)
      character*15 ctstart
      character*80 xlabel

      logical lplot, lcirc, ltrack, lstat, llabel, ltime, datearg
      logical lwedge, lcoast, locean
      real rlon, rlat
      integer symbol, lnblnk
      real*8 pi, raddeg, elev
      parameter (pi = 3.141592654)
      parameter (raddeg = 180.0D0/pi)

      real lon0, lon1, lat0, lat1, latmap0, latmap1, lonmap0, lonmap1

      logical asc, des, ascending

* statinfo variables 

      real*8    eccentricity(3)
      character text*32,sitename*32, occode*6, code*4
      integer*4 noise, wavelength, plate, statinfo

* getorb variables

      integer*4 getorb, status
      real*8 olat, olon, oheight, sec

* odrinfo variables

      integer*4 odrinfo, fd, rep, arc, ver, nrec,
     |          odrt0, odrt1, tstep, begin, end
      real*8    rev
      character satel*8
*

      real hue, lightness, l0, l1, saturation
      real red, green, blue
      integer numcol
      real trackmin, trackmax

      integer projection
      real mapscale, parm1, parm2
      character*10 units

* Initialize

      type='s'
      sat = ""
      cyc = ""
      staposfile='/home/geodyn/models/stations/itrf2000/itrf2000.xyz'
      lresfile='lres.out'
      odrfile='ODR.2'
      device = '/XS'
      alt  = 800.0D3
      elev = 5.0D0
      lcirc=.false.
      ltrack=.false.
      lstat=.false.
      ltime=.false.
      lwedge = .false.
      lcoast = .false.
      locean = .false.
      secb1 = 0.0D0
      secb2 = 2019686400.000D0
      jtrack = 0
      jtotal = 0
      nplotted = 0
      dt = 10.0D0
      asc = .true.
      des = .true.
      lon0 = -180.0
      lon1 =  180.0
      lat0 = -90.0
      lat1 =  90.0

      projection = 5
      mapscale = 0.0
      parm1 = 0.0
      parm2 = 0.0

* Read command-line arguments

      do i=1, iargc()
        call getarg(i,arg)
        if (datearg(arg,secb1,secb2,dt)) then
        elseif (arg(1:4).eq."lat=") then
          read(arg(5:),*) lat0,lat1
        elseif (arg(1:4).eq."lon=") then
          read(arg(5:),*) lon0,lon1
        elseif (arg(1:4).eq."sat=") then
          read(arg(5:),*) sat
        elseif (arg(1:4).eq."cyc=") then
          read(arg(5:),*) cyc
        elseif (arg(1:4).eq."-asc") then
          asc = .true.
          des = .false.
        elseif (arg(1:4).eq."-des") then
          asc = .false.
          des = .true.
        elseif (arg(1:4).eq."odr=") then
          odrfile = arg(5:)
        elseif (arg(1:5).eq."lres=") then
          lresfile = arg(6:)
        elseif (arg(1:4).eq."alt=") then
          read(arg(5:),*) alt
          alt = alt * 1d3
        elseif (arg(1:5).eq."elev=") then
          read(arg(6:),*) elev
        elseif (arg(1:4).eq."dev=") then
          device = arg(5:)
        elseif (arg(1:5).eq."type=") then
          type = arg(6:)
          if (type.eq.'d') then
            trackmin = 0.2
            trackmax = 0.8
          elseif (type.eq.'s') then
            trackmin = 0.0
            trackmax = 8.0
          end if
        elseif (arg(1:2).eq."-s") then
          lstat = .true.
        elseif (arg(1:5).eq."-circ") then
          lcirc = .true.
        elseif (arg(1:2).eq."-g") then
          ltrack = .true.
          lwedge = .true.
        elseif (arg(1:2).eq."-l") then
          llabel = .true.
        elseif (arg(1:2).eq."-t") then
          ltime = .true.
        elseif (arg(1:4).eq."map=") then
          read(arg(5:),*) projection, mapscale, parm1, parm2
        elseif (arg(1:6).eq."-coast") then
          lcoast = .true.
        else
          goto 1300
        endif
      enddo

* Initialize plot and set colors
 
      call pgopen(device)
      call pmdef(projection,mapscale,parm1,parm2)
      call pgsvp (0.10, 0.90, 0.15, 0.95)
      call pmswin(lon0,lon1,lat0,lat1)

      call pgscr(0, 1.0, 1.0, 1.0) ! Background color
      call pgscr(2, 0.7, 0.8, 0.7) ! Continents color
      call pgscr(3, 1.0, 0.0, 0.0) ! Stations color
      call pgscr(15, 0.85, 0.95, 1.0) ! Oceans
*      call pgscr( 3, 1.0, 1.0, 1.0) ! zero stations

      hue = 0.5
      l0 = 0.3
      l1 = 0.8
      saturation = 0.8
      numcol = 100
      lightness = 0.5
      do i=1,numcol
        hue = ((1.0-real(i)/numcol)) * 300.0 + 60.0
        lightness = (1.0-real(i)/numcol) * (l1-l0) + l0
        call grxrgb(hue,lightness,saturation,red,green,blue)
        call pgscr(i+16, red, green, blue)
      end do

      latmap0 = lat0
      latmap1 = lat1
      lonmap0 = lon0
      lonmap1 = lon1
      call pmconv(1,lonmap0,latmap0)
      call pmconv(1,lonmap1,latmap1)
      call pgsci(15)
*      call pgrect(lonmap0,lonmap1,latmap0,latmap1)
*      call pgrect(lon0,lon1,lat0,lat1)
      call pgrect(-1000,1000,-1000,1000)
      call pgsci(2)
      call pmwdb('0.lnd',2,0)
      call pgscir(17,numcol+16)
      call pgsci(1)
      if(type.eq.'d') then
        units='mm/s'
      else
        units='cm'
      end if
      if (lwedge) then
        call pgwedg('BI',1.2,3.0,trackmin,trackmax,units)
      end if

* get station positions and pass start and end times

      write(*,*) 'Loading station coordinate file...'
      call loadcoord (staposfile, mjdref, nstapos, statid, pos, sig)
      write(*,*) 'Reading tracking pass list from '
     | //lresfile(1:lnblnk(lresfile))//'...'
      call readlresout(lresfile,maxpas,station,mean,rms,sec1,sec2,npas)

* get ODR begin and end time

*      status = odrinfo(fd,odrfile,satel,rep,arc,ver,nrec,
*     |             odrt0,odrt1,tstep,begin,end,rev)
*      odrfile = "+" // odrfile
      odrt0 = secb1
      odrt1 = secb2

      if(ltime) then
        call pgsvp (0.1, 0.95, 0.10, 0.3)
        call pgsci(1)
        call pgswin (0.0,real(odrt1-odrt0)/86400.0,0.0,5.0)
        call pgbox ('BCTN', 0.0, 0, 'BCNST', 0.0, 0)
        call chrdat(odrt0,ctstart)
        xlabel = 'days since ' // ctstart
        call pglab(xlabel,'stations','')
      end if

* Build array of tracking start times

      if (.not.ltrack) goto 40
      write(*,*) 'Building array of passes ...'

      sec = max(odrt0,nint(secb1)) * 1.0D0

      nplot = 1
      splot(1) = sec
      nstaplot(1) = 0

  30  continue

      nsta = 0
      ndor = 0
      npra = 0
      nslr = 0

      do i=1,npas
        if(sec.ge.sec1(i) .and. sec.le.sec2(i)) then
          nsta = nsta + 1
          if(station(i).ge.4000.and.station(i).lt.5000) then
            ndor = int(16 + ((rms(i) - trackmin)/trackmax)*numcol)
            if (ndor<=17) ndor = 17
            if (ndor>=116) ndor = 116 
!            if(rms(i).le.trackmin) ndor = 17
!            if(rms(i).ge.trackmax) ndor = 17 + numcol
          elseif(station(i).ge.7700.and.station(i).lt.7799) then
            npra = npra + 1
          else
            nslr = int(17 + ((abs(rms(i)) - trackmin)/trackmax)*numcol)
            if(rms(i).le.trackmin) nslr = 17
            if(rms(i).ge.trackmax) nslr = 16 + numcol
          endif
        endif
      enddo

      if (type.eq.'d') then 
        i=ndor
      elseif (type.eq.'s') then 
        i=nslr
      elseif (type.eq.'pd') then
        i=npdop
      elseif (type.eq.'pr') then
        i=npran
      elseif (type.eq.'p') then
        i=npra
      elseif (type.eq.'a') then
        i=nsta
      else
        i=0
      endif

      if (i.ne.nstaplot(nplot)) then
        nplot = nplot + 1
        nstaplot(nplot) = i
        splot(nplot) = sec
      endif

      if(i.gt.0) jtrack = jtrack + 1

      jtotal = jtotal + 1

      if(ltime) then
        call pgpoint(1,real(sec-odrt0)/86400.0,real(i),1)
      endif

      sec = sec + dt

      if(sec.gt.odrt1*1.0D0.or.sec.gt.secb2) goto 35

      goto 30

 35   continue

      nplot = nplot + 1
      nstaplot(nplot) = 0
      splot(nplot) = sec
      
      write(*,*) 'Plotting passes on map ...'

      do k=1, 255      
      do i=1,nplot-1
      if(nstaplot(i).eq.k) then
        npoints = int((splot(i+1)-splot(i))/dt)
        do j=1, npoints
          status = getorb(splot(i)+j*dt,
     |              olat,olon,oheight,odrfile,.true.)
          x(j) = real(olon)
          y(j) = real(olat)
        enddo

        ascending = .false.
        if(y(1).gt.y(npoints)) ascending = .true.

        if((ascending.and.asc).or.(.not.ascending.and.des)) then
          call pgsci(nstaplot(i))
          if(nstaplot(i).ge.1) then
            call pgslw(5)
          else
            call pgslw(2)
          endif
          call pmconv(npoints,x,y)
          call mapline(npoints,x,y)
        endif
      endif

      enddo
      enddo

      write(*,*) 
     | (real(jtrack) / real(jtotal)) * 100E0, '% tracking coverage.' 

 40   continue

      call pgslw(1)

      if (.not.(lcirc.or.lstat)) goto 50

* plot range circles for stations

      do j=1,npas
        if (nplotted.eq.0) then
          nplotted = 1
          splotted(nplotted) = station(j)
        else
          do k=1,nplotted
            if (splotted(k).eq.station(j)) goto 45
          enddo
        end if
        nplotted = nplotted + 1
        splotted(nplotted) = station(j)
        do i=1,nstapos
        lplot = .false.
        if (statid(i).eq.station(j)) then
          if (station(j).ge.4000.and.station(j).lt.5000) then
            if (type.eq.'d') lplot = .true.
            if (type.eq.'a') lplot = .true.
          elseif (station(j).lt.7799.and.station(j).ge.7700) then
            if (type(1:1).eq.'p') lplot = .true.
            if (type.eq.'a') lplot = .true.
          else
            if (type.eq.'s') lplot = .true.
            if (type.eq.'a') lplot = .true.
          end if
          if (lplot) then 
           mjd = sec1(j) / 86400.0D0 + 46066.0D0
*           print *, j, station(j), sitename, occode, code
           if(statinfo(mjd,station(j),sitename,occode,code,
     |                   noise,wavelength,plate,eccentricity).eq.1) 
     |     then
              print *, 'pass, station: ', j, station(j)
              stop 'station not found'
           endif
           call pgsci(1)
           symbol = -3
           xyz(1)=pos(1,i)
           xyz(2)=pos(2,i)
           xyz(3)=pos(3,i)
           call xyzgeo (xyz, r, lat, lon, height)
*           write(*,'(i4,1x,3f20.4)') statid(i), lat, lon, height
           rlon = real(lon*raddeg)
           rlat = real(lat*raddeg)

* Plot station marker

           if (lstat) then 
            if (rlon > lon0 .and. rlon < lon1 .
     &     and. rlat > lat0 .and. rlat < lat1) then
             call pgsch(1.50)
             call pgsci(3)
             call pmconv(1,rlon,rlat)
             call pgpoint(1,rlon,rlat,symbol)
             call pgsci(1)
             if (llabel) then
               text = occode
               call pgsch(0.8)
               call pgtext(rlon+1.0,rlat+1.0,text)
             endif
            endif
           endif
           if (lcirc) then
            call trackcir(alt,elev,lon*raddeg,lat*raddeg,height*1d-3)
           endif
         endif
        endif
        enddo
 45     continue
      enddo

 50   continue

      call pgsci(1)
      if (lcoast) call pmwdb('1.cil',1,1)
      IF (projection .ge. 40) THEN 
        call pmbox('H-',720.0,0,'H-',360.0,0)
      ELSE
        call pmbox('IBC-TS',30.0,2,'IBC-TS',30.0,2)
      END IF

      IF (type=='d' .AND. LEN(sat)>0 .AND. LEN(cyc)>0) THEN
        title = "DORIS residuals for " // sat // " cycle: " // cyc
        call pgsch(1.50)
        CALL pglab("","",title)
        call pgsch(1.00)
      END IF
      call pgend

      goto 9999

1300  write (0,600)
600   format('trackingmap - makes plot of tracking data
     |  using ODR and lres.out files'/
     |'syntax: trackingmap [ options ]'//
     |'where [ options ] are:'/
     |'t=t0,t1     : Specify time interval (MJD,YYMMDD,SEC85)',
     |' and step size (sec, def: 60)'/
     |'              ... or use mjd=, doy=, ymd=, sec='/
     |'-s          : Plot station locations'/
     |'-l          : Plot text labels next to stations (-s reqd)'/
     |'-c          : Plot visibility circles'/
     |'-g          : Plot groundtracks'/
     |'-t          : Plot timeline'/
     |'odr=name    : Specify odr filename (default: ODR.2)'/
     |'lres=name   : Specify lres output file (default: lres.out)'/
     |'type=x      : Specify ploting of tracking type x, where x is'/
     |'              s (SLR), d (DORIS), p (PRARE), a (all)')
9999  end

      subroutine readlresout(filenm,maxpas,station,mean,rms,
     |                       tindx1,tindx2,n)

      real*8 date, hour, minute, second, time1, time2
      real*8 datetime1, datetime2
      real*8 sec1, sec2, sec85, tindx1(*), tindx2(*), mean(*), rms(*)
      integer*4 station(*), n, unit, maxpas
      character*80 filenm, line

      n=0

      if(filenm(1:1).eq.' ') then
        unit=-1
      else
        unit=10
        open(unit,file=filenm)
      end if

 10   continue
      read(unit,'(a)',end=20) line
      if (n.gt.maxpas) then
        write(*,*) 'Maximum number of passes reached in lresread'
        goto 20
      end if
      if (line(1:6).eq."      ") goto 10

      n = n + 1
      read(line(37:40),*) station(n)
      read(line(63:69),*) mean(n)
      read(line(71:77),*) rms(n)

      read(line(1:6),*) date
      read(line(8:9),*) hour
      read(line(11:12),*) minute
      read(line(14:15),*) second
      time1 = hour*1d4+minute*1d2+second
      read(line(19:20),*) hour
      read(line(22:23),*) minute
      read(line(25:26),*) second
      time2 = hour*1d4+minute*1d2+second

      datetime1 = date*1d6+time1
      datetime2 = date*1d6+time2
      sec1 = sec85(4,datetime1)
      sec2 = sec85(4,datetime2)

      if (time2.lt.time1) sec2=sec2+86400.0D0
      tindx1(n) = sec1
      tindx2(n) = sec2

      goto 10

  20  continue
      close(unit)
      end

      subroutine trackcir(alt,elev,lon,lat,r)

      real*8 alt, elev, lon, lat
      real*8 pi, ae, finv, alpha, gamma, x
      real*8 r
      real rlon,rlat,radius
      integer np, k
      parameter (np=100)
      real cx(np), cy(np), lx(2), ly(2)
      parameter (pi = 3.141592654)

      call getearth(ae, finv)

      rlon=real(lon)
      rlat=real(lat)
      r = r + ae

      gamma = pi / 2 + elev * (pi/180)
      x=0.5*(-2*r*cos(gamma)-sqrt((2*r*cos(gamma))**2
     |             +8*r*alt+4*alt**2))
      alpha = acos(-(x**2-r**2-(r+alt)**2)/(2*r*(r+alt)))

      radius = real(alpha * ae * 1e-3)
      call rngcir(rlon,rlat,radius,np,0.0,359.9,cx,cy)
      call pmconv(np,cx,cy)
      call pgsci(1)
      do k=1,np-1
        lx(1)=cx(k)
        lx(2)=cx(k+1)
        ly(1)=cy(k)
        ly(2)=cy(k+1)
        if(lx(1).lt.-180.0.or.lx(2).lt.-180.0) then
          lx(1) = lx(1) + 360.0
          lx(2) = lx(2) + 360.0
          call pgline(2,lx,ly)
        elseif(lx(1).ge.180.0.or.lx(2).ge.180.0) then
          lx(1) = lx(1) - 360.0
          lx(2) = lx(2) - 360.0
          call pgline(2,lx,ly)
        else
          call pgline(2,lx,ly)
        endif
      enddo
C      call pgline(np,cx,cy)

      end

      subroutine mapline(npoints,x,y)

      integer*4 npoints
      real x(npoints), y(npoints),xx(2),yy(2)

      integer*4 i

      do i=1,npoints-1
        if(abs(x(i)-x(i+1)).gt.100.0) then
        else
          xx(1)=x(i)
          xx(2)=x(i+1)
          yy(1)=y(i)
          yy(2)=y(i+1)
          call pgline(2,xx,yy)
        endif
      enddo

      end
