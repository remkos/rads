      subroutine update
      implicit none
      include "data.inc"
      include "stat.inc"
      integer it,icol,ip,vpnt

      do it=1,ngood

         do icol=1,npar
            ip=vpnt(it,icol)
            param(ip)=param(ip)+vector(ip)
            adjst(ip)=vector(ip)
         enddo

      enddo

      end
