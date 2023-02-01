      program test_perth

*  Test PERTH

      implicit double precision (a-h,o-z)
      logical select(8), isdata

      write(6,17)
      time = 49046.803493 d0
      dlat = -21.00
      dlon = 5.0
      do i=1,5
         call perth3( dlat,dlon,time,tide,isdata )
         if (isdata) then
            write(6,15) dlat,dlon,time,tide,isdata
         else
            write(6,16) dlat,dlon,time,isdata
         endif
         dlon = dlon + 5.d0
      end do
      dlon = 5.0
      do i=1,8
         call perth3( dlat,dlon,time,tide,isdata )
         write(6,15) dlat,dlon,time,tide,isdata
         time = time + 1.d0/24.d0
      end do
 15   format(1x,2f10.2,f15.6,f10.2,6x,l1)
 16   format(1x,2f10.2,f15.6,16x,l1)
 17   format(5x,'latitude  longitude',7x,'time',6x,'tide    isdata?'/)
      stop
      end
