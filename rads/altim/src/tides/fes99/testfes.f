c########################################################################
c      
      program testfes
c
c------------------------------------------------------------------------
c  ROUTINE : testfes.f
c
c  DESCRIPTION : test the tide prediction algorithm
c
c  PROGRAMMERS : F. LEFEVRE
c 
c  DATE : 09/11/2000
c------------------------------------------------------------------------
        integer*4      jnasa, i, istat, iasc
        real*4         lon, lat
        real*4         h, hlp, hload, tide, tideload
        real*8         heure
        character*20   model
        character*512  path

        model='fes99'
        path='../data/'

        lon = -7.688
        lat = 59.195

        jnasa=9131
        heure=0.d0

c--- Use binary files : iasc=1
        iasc=0
        
        do i=0,24
          call fes_tide(model,path,h,hlp,hload,lat,lon,
     &     jnasa,heure,iasc,istat)
          if(i.eq.0) then
            write(*,10) 'NASA_d','Hour','Latitude','Longitude',
     &  'Short_tid','LP_tid','Pure_Tide','Geo_Tide','Rad_Tide'
          endif
          tide=h+hlp
          tideload=h+hlp+hload
          write(*,20) jnasa,heure,lon,lat,h,hlp,tide,tideload,hload
          heure=heure+1.d0
        enddo

   10 format(a6,1x,a5,7(1x,a9))
   20 format(i6,1x,f5.2,7(1x,f9.3))
   
      end
