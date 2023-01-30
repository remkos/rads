c
c     FUNCTION RLAM computes longitude for given declination and tracknr.
c
      function rlam(it,decl)
      include "maxcom.inc"
      integer*4 it
      real*8 sinl,decl,rlam,u,arglat
      sinl=tan(decl)/tan(inc(it))
      if (sinl.gt.+1) then
         sinl=+1
      else if (sinl.lt.-1) then
         sinl=-1
      endif
      rlam=asin(sinl)
      u=arglat(decl,inc(it))
      rlam=node(it)+rlam-rotate/theta(it)*u
      end
