* Test driver for perth2
      double precision dlat,dlon,time,tide,lptide
      logical          isdata
*
      write(6,12)
      dlat = -42.
      dlon = 330.
      time = 49100.
      call perth2( dlat,dlon,time,tide,isdata )
      call lpeqmt( time*86400.d0, dlat, lptide )
      tide = tide + lptide
      write(6,11) time,dlat,dlon,tide,isdata
*
      time = 49100.5
      call perth2( dlat,dlon,time,tide,isdata )
      call lpeqmt( time*86400.d0, dlat, lptide )
      tide = tide + lptide
      write(6,11) time,dlat,dlon,tide,isdata
 11   format(F12.5,2f8.2,f8.2,6x,L1)
 12   format(5x,'time',7x,'lat',5x,'lon',5x,'tide')
      stop
      end
