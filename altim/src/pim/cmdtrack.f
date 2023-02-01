      subroutine cmdtrack
      include "pim.inc"
      real lat0,lat1,lon0,lon1,sdt
      real*8 sec85,dt0,dt1,dt,dumd
      integer it0,it1,asc,unit,freeunit
      logical dpimopt,pimopt,l,chropt
      character*80 filenm

* Draw groundtracks

      call pgsave
10    if (pop1cmd('TRACK',argum)) then
	 ci=1
	 ls=1
	 lw=1
	 asc=0
	 dt=60d0
	 lon0=xw0
	 lon1=xw1
	 lat0=yw0
	 lat1=yw1
	 name='ers1'
	 if (pimopt('-asc',argum,dum,dum,dum,dum)) asc=1
	 if (pimopt('-des',argum,dum,dum,dum,dum)) asc=-1
	 l=pimopt('ci=',argum,ci,dum,dum,dum)
	 l=pimopt('ls=',argum,ls,dum,dum,dum)
	 l=pimopt('lw=',argum,lw,dum,dum,dum)
         if (pimopt('dt=',argum,sdt,dum,dum,dum)) dt=sdt
         l=pimopt('lat=',argum,lat0,lat1,dum,dum)
         l=pimopt('lon=',argum,lon0,lon1,dum,dum)
	 l=chropt('sat=',argum,name)
	 if (name.eq.'geosat') then
	    dt0=870617154520.0d0
	    dt1=870704165720.0d0
	 else if (name.eq.'ers1' .or. name.eq.'ers1:35') then
	    name='ers1'
            dt0=930101d0
            dt1=930205d0
	 else if (name.eq.'ers1:3') then
	    name='ers1'
            dt0=940101d0
            dt1=940104d0
	 else if (name.eq.'ers1:168') then
	    name='ers1'
            dt0=940410d0
            dt1=940928d0
	 else if (name.eq.'ers1:168b') then
	    name='ers1'
            dt0=940928d0
            dt1=950321d0
	 else if (name.eq.'topex') then
            dt0=930130.054873d0
            dt1=dt0+9.916d0
	 else if (name.eq.'ers2') then
	    dt0=950601d0
	    dt1=950706d0
	 endif
	 call pgsls(nint(ls))
	 call pgslw(nint(lw))
	 call pgsci(nint(ci))
	 write (0,'("Plotting groundtracks ...")')
	 if (chropt('file=',argum,filenm)) then
	    unit=freeunit()
	    open (unit,file=filenm,status='old')
15	    read (unit,*,end=20) dt0,dt1
            it0=nint(sec85(0,dt0))
            it1=nint(sec85(0,dt1))
	    call tracks(name,lat0,lat1,lon0,lon1,it0,it1,nint(dt),asc)
	    goto 15
20	    continue
	 else
            l=dpimopt('t=',argum,dt0,dt1,dt,dumd)
            it0=nint(sec85(0,dt0))
            it1=nint(sec85(0,dt1))
	    call tracks(name,lat0,lat1,lon0,lon1,it0,it1,nint(dt),asc)
	 endif
	 goto 10
      endif
      call pgunsa
      end

      subroutine tracks(sat,lat0,lat1,lon0,lon1,t0,t1,dt,asc)
      include "pim.inc"
      character*(*) sat
      character*80 filenm
      real lat0,lat1,lon0,lon1,lato,lono
      real*8 lat,lon,orb
      integer t,t0,t1,dt
      integer l,lnblnk,getorb,asc
      logical move

      filenm='/user/altim'
      call checkenv('ALTIM',filenm,l)
      if (sat.eq.'geosat') then
	 filenm(l+1:)='/data/ODR.GEOSAT'
      else if (sat.eq.'topex') then
	 filenm(l+1:)='/data/ODR.TOPEX/JGM-2'
      else if (sat.eq.'ers1') then
	 filenm(l+1:)='/data/ODR.ERS-1/dgm-e04'
      else if (sat.eq.'ers2') then
	 filenm(l+1:)='/data/ODR.ERS-2/latest'
      else
         STOP "PIM: unknown satellite in CMDTRACK"
      endif

      move=.true.
      do t=t0,t1,dt
	 lato=lat
	 lono=lon
	 l=getorb(dble(t),lat,lon,orb,filenm,.true.)
         if ((lat.le.lat0 .or. lat.ge.lat1 .or.
     |       lon.le.lon0 .or. lon.ge.lon1) .and.
     |       (lato.le.lat0 .or. lato.ge.lat1 .or.
     |       lono.le.lon0 .or. lono.ge.lon1)) move=.true.
	 if (lat.gt.lato .and. asc.lt.0) move=.true.
	 if (lat.lt.lato .and. asc.gt.0) move=.true.
	 if (l.gt.0) stop 'pim: error using getorb'
	 if (lon.lt.lon0) lon=lon+360
	 if (lon.gt.lon1) lon=lon-360
	 if (move) then
	 else if (project.gt.40) then
	    if (lat.lt.lat0 .or. lat.gt.lat1) move=.true.
	 else if (lon-lono.gt.180) then
	    x=lon-360
	    y=lat
	    call pmconv(1,x,y)
	    call pgdraw(x,y)
	    x=lono+360
	    y=lato
	    call pmconv(1,x,y)
	    call pgmove(x,y)
	 else if (lono-lon.gt.180) then
	    x=lon+360
	    y=lat
	    call pmconv(1,x,y)
	    call pgdraw(x,y)
	    x=lono-360
	    y=lato
	    call pmconv(1,x,y)
	    call pgmove(x,y)
         endif
	 x=lon
	 y=lat
	 call pmconv(1,x,y)
	 if (move) then
	    call pgmove(x,y)
	 else
	    call pgdraw(x,y)
	 endif
	 move=.false.
      enddo

      end
