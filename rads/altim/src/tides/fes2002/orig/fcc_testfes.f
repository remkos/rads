c########################################################################
c      
      program fcc_testfes
c
c------------------------------------------------------------------------
c  PROGRAM : fcc_testfes.f
c
c  DESCRIPTION : test the tide prediction algorithm for FES2002
c                fortran F77 version
c
c  PROGRAMMERS : F. LEFEVRE, F. BRIOL
c 
c  DATE : June, 3rd 2003
c------------------------------------------------------------------------
c  rc      : return code (problem if rc != 0)
c  lat     : latitude
c  lon     : longitude
c  time    : time in CNES Julian days
c  hour    : hour
c  tide    : short tides (semi_diurnal and diurnal tides)
c  load    : loading effects
c  lp      : long period tides
c
c  tide+lp         = pure tide (as seen by a tide gauge)
c  tide+lp+load    = geocentric tide ((as seen by a satellite)
c  CNES Julian day = NASA Julian day + 2922
c  CNES Julian day 0 is at midnight between the 31 December 
c                          and 01 January 1950 AD Gregorian
c
c########################################################################

      integer*4 rc
      real*8	lat
      real*8	lon
      real*8	time
      real*8	tide
      real*8	load
      real*8	lp
      
      integer*4	hour
      integer*4	fes2002
      character*512	dir
      
      dir="../data"
      fes2002 = 1
      lat	= 59.195d0
      lon	= -7.688d0
c    Time : 1983-01-01 00:00:00.0
      time	= 12053d0

c    Initialize memory for FES algorithms
      call fcc_initfes(dir, fes2002, rc)
      if(rc .eq. 1) stop
      
      write(*,10) 'JulDay','Hour','Latitude','Longitude',
     &  'Short_tid','LP_tid','Pure_Tide','Geo_Tide','Rad_Tide'

      do hour=0,23

c      Compute tide
        call fcc_festide(lat, lon, time, tide, load, lp, rc)
        if(rc .eq. 1) stop
	

        write(*,20)	time,
     & 			hour, 
     & 			lon, 
     & 			lat,
     & 			tide,
     & 			lp,
     & 			tide+lp,
     & 			tide+lp+load,
     & 			load
	
	time = time + 1 / 24d0

      enddo

c    Free memory for FES algorithms
      call fcc_freefes()

   10 format(a12,1x,a5,7(1x,a9))
   20 format(f12.5,1x,i5,7(1x,f9.3))

      end

