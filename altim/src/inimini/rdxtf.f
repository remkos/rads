* rdxtf now takes care of ALL checks on good(it) !
* In all cases itr(it) = 0  where  it is the normal track-number (as
* stored in the file).  itr is the track pointer (sequential in time,
* continuous, and only good tracks).
*
* Note
* it=itr(j)     gives pointer it for tracknr j
* j=itrinv(it)  gives tracknr j for pointer it
*
      subroutine rdxtf
      implicit none
      include "data.inc"
      include "satcat.inc"
      include "stat.inc"
      include "init.inc"
      integer*4	i,ip,icol,irec,l,iflags,it,vpnt,inc,sal,epasc,lonasc,
     |		startt,stopt,ipar(maxpar),isig(maxpar)
      integer*2 jobs,satid,nxor,noalt,flags
      real*8	toff(10)/0d0,0d0,3197d0,7*0d0/

* Open track catalogue

      open (20,file=xtf,form='unformatted',access='direct',
     |         recl=xtflen,status='old',err=1000)

* Start reading catalog, depends on filetype

      ngood=0
      do irec=2,ntrack+1
         read (20,rec=irec) jobs,satid,nxor,noalt,inc,sal,
     |                      epasc,lonasc,startt,stopt,
     |                      (ipar(i),i=1,npar),
     |                      (isig(i),i=1,npar),flags
         it=wrap(int(jobs))
         itr(it)=0
         t(it)=startt+toff(satid)*86400
	 u0(it)=sal

* Set flags and pointers

         iflags=flags

* The pass is good if it is marked good in the catalog,
* and if the satellite index was requested and if there
* are at least 2*NPAR residuals

	 if (iand(iflags,512).ne.0 .and. inimode.le.0) then
	 else if (iand(iflags,256).eq.0) then
	 else if (inimode.eq.1 .and. nxor.lt.2*npar) then
	 else if (inimode.eq.2 .and. noalt.lt.2*npar) then
	 else
            ngood=ngood+1
            itrinv(ngood)=it

* Store satellite ID temporarily in t_nres

	    t_nres(it)=satid
            do icol=1,npar

* A priori parameters are stored in vector
* A priori sigmas are stored in matrix
* ... but just until all tracks have been sorted

               ip=vpnt(it,icol)
	       if (.not.init) then
	          vector(ip)=ipar(icol)/1e6
	       else if (icol.eq.1) then
		  vector(ip)=-abias(satid)
	       else
		  vector(ip)=0e0
	       endif
            enddo
            matrix(it)=isig(1)/1e6
	    if (init) matrix(it)=orberr(satid)
         endif
      enddo

      return

 1000 l=index(xtf,' ')-1
      write (6,2000) xtf(1:l)
      call fin('inimini: stop')

 2000 format ('Track catalog does not exist and I really ',
     |        'absolutely need it'/'Filename : ',a)
      end
