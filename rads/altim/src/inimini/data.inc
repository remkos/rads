* data.inc

* maxpar        maximum number of parameters per track
* mtrack        maximum nuber of tracks
* maxsiz        maximum size of the normal matrix * 8 bytes
* bndwth	width of the bandmatrix
* bndlng	length of the bandmatrix = npar*ngood
* itr		conversion data track number -> internal track number
* itrinv	conversion internal track number -> data track number
* t_nrestot	total number of residuals per track
* t_nres	number of used residuals per track
* t_rres	rms of residuals per track
* t_eres	weighted rms of residuals per track
* t		time of the track
* u0		start argument of latitude of the track (microdeg)
* good		quality per track
* satel		satellite id per track
* ngood		number of apriori valid tracks
* param		orbit error parameters
* sigma		orbit error sigmas
* adjst		latest parameter adjustment

      integer maxpar,mtrack,maxsiz
      parameter (maxpar=5,mtrack=65536,maxsiz=5000000)

      real    param(mtrack*maxpar),adjst(mtrack*maxpar),sigma(mtrack),
     |        t(mtrack)
      integer itr(mtrack),itrinv(mtrack),u0(mtrack),t_nrestot(mtrack),
     |        t_nres(mtrack),satel(mtrack),ngood,
     |	      bndwth,bndlng
      real    t_rres(mtrack),t_eres(mtrack)
      logical good(0:mtrack)

      common /tracks/ t,itr,itrinv,good,u0,t_nrestot,t_nres,t_eres,
     |		t_rres,satel,param,sigma,adjst,ngood,
     |	        bndwth,bndlng

* matrix	band matrix of normal equations
* vector	right-hand side vector

      real    matrix(maxsiz),vector(mtrack*maxpar)
      common /equations/ matrix,vector

* xxf,xtf,xaf	XXF and XTF and XAF file names
* xxflen	length of xover data record
* xtflen	length of track data record
* extend	extended data records
* short		short data records
* ntrack	number of tracks in XTF
* nrec		number of xovers in XXF
* npar		number of parameters per track
* init		initialise parameters at first iteration
* inimode	mode: 1=inimini, 2=inisurf, 0=iniedit

      character	xxf*80,xtf*80,xaf*80
      integer   xxflen,xtflen,xaflen,ntrack,nrec,npar,inimode
      logical   extend,short,init

      common /files/ xxf,xtf,xaf,xxflen,xtflen,xaflen,
     |		ntrack,nrec,extend,short,init,inimode,npar
