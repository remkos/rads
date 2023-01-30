c########################################################################
c      
      program fes
c
c------------------------------------------------------------------------
c ROUTINE : fes.f
c
c DESCRIPTION : tide prediction algorithm for FES99
c
c PROGRAMMERS : F. LEFEVRE
c 
c DATE : 21/11/2000
c------------------------------------------------------------------------
c
c Variables :
c -----------
c lon      = longitude
c lat      = latitude
c jnasa    = date is given as on AVISO/CDROM
c h        = pure diurnal and semi diurnal tide elevation at a given time
c hlp      = pure long period tide elevation at a given time
c hload    = radial loading tide at a given time
c tide     = pure tide elevation at given time (tide=h+hlp)
c tideload = geocentric tide elevation at given time (tideload=h+hlp+hload)
c iasc     = 1 for ascii files
c          = 0 for binary files
c
c Subroutines :
c -------------
c convert_ymd2nasa : convert a date yy/mm/dd to NASA day to be applied in
c                    our algorithm
c------------------------------------------------------------------------
        integer*4      jnasa, i, istat, iasc
        integer*4      iyear, imonth, iday
	real*4         lon, lat
	real*4         h, hlp, hload, tide, tideload
	real*8         hour
	character*20   model
	character*512  path

        model='fes99'
	path='../data/'

c-- Coordinates
	lon = 303.419527
	lat = 63.983873

c--- Use binary files : iasc=0
        iasc=0

c-- NASA day is used in fes_tide subroutine to compute tide elevation
c Date is given in "Jnasa" (as on AVISO/CDROM) and the time is given
c in hours, during the day (double precision, between 0.d0 and 24.d0)
c For instance Jnasa=12692 for 1 October ,1992.

c-- Convert NASA day in dd/mm/yyyy
c	jnasa=9131
c       call calend1(jnasa+2922,iday,imonth,iyear)
c       print*, iday, imonth, iyear


c-- Enter date yy/mm/dd
c   For instance : iyear=1983 imonth=1 iday=1 <=> January 1st, 1983 00:00
        iyear  = 2000
        imonth = 3
        iday   = 23
	hour  = 15.0669444d0

c-- Convert date yy/mm/dd to NASA day
c This subroutine allows to give useful date format
c which is converted in NASA day
        call convert_ymd2nasa(jnasa,iyear,imonth,iday)
	print*, 'jnasa = ', jnasa

        do i=0,0
	  call fes_tide(model,path,h,hlp,hload,lat,lon,
     &                  jnasa,hour,iasc,istat)
          if(i.eq.0) then
	    write(*,10) 'Year,','Mo','Da','Hour','Latitude','Longitude',
     &  'Short_tid','LP_tid','Pure_Tide','Geo_Tide','Rad_Tide'
          endif
          tide=h+hlp
          tideload=h+hlp+hload
          call calend1(jnasa+2922,iday,imonth,iyear)
	  write(*,20) iyear,imonth,iday,hour,lon,lat,
     &                h,hlp,tide,tideload,hload
	  hour=hour+1.d0
          if (hour.ge.24d0) then
             hour=hour-24d0
             jnasa=jnasa+1
          endif
	enddo

   10 format(a4,1x,a2,1x,a2,1x,a5,7(1x,a9))
   20 format(i4,1x,i2,1x,i2,1x,f5.2,7(1x,f9.3))
	
        end
