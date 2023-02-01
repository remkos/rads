*
* Compute argument of latitude from declination and inclination
*
      function arglat(decl,incl)
      real*8 arglat,decl,incl,sinu
      sinu=sin(decl)/sin(incl)
      if (sinu.gt.+1.0d0) then
         sinu=+1.0d0
      else if (sinu.lt.-1.0d0) then
         sinu=-1.0d0
      endif
      arglat=asin(sinu)
      end
