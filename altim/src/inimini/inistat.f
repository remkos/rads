      subroutine initstat
      implicit none
      include "stat.inc"
      include "data.inc"
      integer ix

* Initialise statistics per iteration


* Initialise xover rms per satellite combination

      do ix=1,maxcmb
         s_rres(ix)=0e0
         s_nres(ix)=0
      enddo

* Initialise total correction mean and rms per satellite

      h(1)=0
      do ix=1,maxsat
         s_mcor(ix)=0e0
         s_rcor(ix)=0e0
	 s_adjust(ix)=0e0
         s_ncor(ix)=0
	 h(ix)=h(ix-1)+(ix-1)
      enddo

* Initialise number of xovers per track

      do ix=1,ngood
         t_nres(ix)=0
         t_nrestot(ix)=0
         t_rres(ix)=0e0
         t_eres(ix)=0e0
      enddo

      end
