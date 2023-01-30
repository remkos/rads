* file: piggy.inc
      integer mx,nx,ny
      parameter (mx=1201*1201)
      real*8 h(mx),lon0,lon1,dlon,lat0,lat1,dlat
      common /cpiggy/ h,lon0,lon1,dlon,lat0,lat1,dlat,nx,ny
