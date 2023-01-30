      function mpnt(track1,track2,par1,par2)
      implicit none
      integer track1,track2,par1,par2,ip1,ip2,p1,p2,mpnt,vpnt
      include "data.inc"

* Determine matrix pointer mpnt.
* Return mpnt=0 when "out of band".

      if (abs(track1-track2).ge.bndwth/npar) then
	 mpnt=0
      else
	 ip1=vpnt(track1,par1)
	 ip2=vpnt(track2,par2)
	 p1=min(ip1,ip2)
	 p2=max(ip1,ip2)

* Band storage mode used by LAPACK`s SPBTRF(`U`)

	 mpnt=p2*(bndwth-1)+p1
      endif
      return
      end
