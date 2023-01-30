* stat.inc

* maxsat        maximum number of satellites
* maxcmb	maximum number of sat/sat combinations
* resrms	total xover rms
* resrms0	total xover rms zerost iteration
* adjust	maximum orbit adjustment
* s_mcor	mean orbit correction per satellite
* s_rcor	rms orbit correction per satellite
* s_adjust	maximum orbit adjustment per satellite
* s_ncor	number of corrections per satellite
* s_nres	number of residuals per sat/sat combination
* s_xormns	rms residuals per sat/sat combination
* iter		number of iterations
* ntrks		number of tracks used
* h		pointer

      integer maxsat,maxcmb
      parameter (maxsat=20,maxcmb=maxsat*(maxsat+1)/2)

      real*8  	s_mcor(maxsat),
     |		s_rcor(maxsat),s_rres(maxcmb)
      real      s_adjust(maxsat),resrms,resrms0,adjust
      integer*4	s_ncor(maxsat),s_nres(maxcmb),nrestot,ntrks,
     |		iter,h(maxcmb)

      common /cstat/ s_mcor,s_rcor,s_rres,s_adjust,
     |	    resrms,resrms0,adjust,iter,nrestot,ntrks,s_ncor,s_nres,h
