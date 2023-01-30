      real*8 lat,lon,h,xyz(3),r
      read (5,*) lat,lon,h
      call geoxyz(lat,lon,h,xyz,r)
      write (*,60) lat,lon,h,r,xyz
60    format (7f12.3)
      end
