*************************************************************
      function wieber(kx,nxc,ky,nyc,pi)
      integer kx,ky,nxc,nyc
      real x,y,wieber,pi
*
* Smoothing function for overlaps
*
      x=real(kx)/nxc
      y=real(ky)/nyc
      if (x+y.eq.1.) then
         wieber=.5
      else if (y.gt.x) then
         wieber=.5*(1+cos(x*pi/(1-y+x)))
      else
         wieber=.5*(1+cos(y*pi/(1-x+y)))
      endif
      end
