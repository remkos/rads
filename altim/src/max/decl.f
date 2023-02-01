c
c     FUNCTION DECL compute declination for given longitude and tracknr.
c
      function decl(it,rlam)
      include "maxcom.inc"
      integer*4 it
      real*8 fact,rlam,xlam,u,cosinc,decl
      fact=rotate/theta(it)
      cosinc=cos(inc(it))
      xlam=rlam-node(it)
      if (xlam.ge.+pi/2) then
         u=-pi/2
      else if (xlam.le.-pi/2) then
         u=+pi/2
      else
         u=atan(tan(xlam)/cosinc)
      endif
      xlam=rlam-node(it)+fact*u
      u=atan(tan(xlam)/cosinc)
      xlam=rlam-node(it)+fact*u
      u=atan(tan(xlam)/cosinc)
      xlam=rlam-node(it)+fact*u
      u=atan(tan(xlam)/cosinc)
      xlam=rlam-node(it)+fact*u
      decl=atan(tan(inc(it))*sin(xlam))
      end
