      function vpnt(track,par)
      implicit none
      include "data.inc"
      integer track,par,vpnt

      vpnt=(track-1)*npar+par
      end
