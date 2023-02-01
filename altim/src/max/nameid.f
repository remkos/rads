c
c
c     INTEGER FUNCTION NAMEID switches from satellite name to satellite index
c
c
      function nameid(satnam)
      integer*4 k,nameid
      character*8 satnam
      include "maxcom.inc"
      do k=1,maxsat
         if (satnam.eq.name(k)) then
            nameid=k
            return
         endif
      enddo
      stop 'max: ADR satellite not in name table'
      end
