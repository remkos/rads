      subroutine sort
      implicit none
      include "data.inc"
      integer*4 i,j,j1,j2,it,ip,icol,vpnt

* First sort the trackpointers in chronological order
* it=itr(j)    gives pointer it for tracknr j
* j=itrinv(it) gives tracknr j for pointer it

      do i=1,ngood
	 do it=1,ngood-i
            j1=itrinv(it)
            j2=itrinv(it+1)
	    if (t(j1).gt.t(j2)) then
               itrinv(it)=j2
               itrinv(it+1)=j1
	    endif
	 enddo
      enddo

      do it=1,ngood
         j=itrinv(it)
         itr(j)=it
	 good(it)=.true.
	 satel(it)=t_nres(j)
      enddo

* Now put all the a priori parameters and sigmas in the right places

      do it=1,ngood
	 j=itrinv(it)
	 do icol=1,npar
	    ip=vpnt(it,icol)
	    i =vpnt( j,icol)
	    param(ip)=vector(i)
	 enddo
	 sigma(it)=matrix(j)
      enddo

      end
