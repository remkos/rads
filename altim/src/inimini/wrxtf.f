      subroutine wrxtf
      implicit none
      include "data.inc"
      include "grid.inc"
      include "stat.inc"
      include "init.inc"
      integer i,irec,ip,icol,it,vpnt

      character	ftype*4
      integer*4 inc,sal,epasc,lonasc,startt,stopt,
     |          ipar(maxpar),isig(maxpar)
      integer*2 jobs,satid,nxor,noalt,flags


* Open track catalog

      open (20,file=xtf,form='unformatted',access='direct',
     |      recl=xtflen,status='old',err=1000)
      read (20,rec=1) ftype

* Start reading catalog

      do irec=2,ntrack+1
         read (20,rec=irec) jobs,satid,nxor,noalt,inc,sal,
     |                      epasc,lonasc,startt,stopt,
     |                      (ipar(i),i=1,npar),
     |                      (isig(i),i=1,npar),flags
         it=itr(wrap(int(jobs)))

* Tracks are marked invalid when:
* - not used in matrix (it=0)
* - number of used crossovers <= 2*npar | Indicated by good(it)
* - edit factor > edttrk                | See subroutine wrstat()

	 if (it.gt.0) then
	    do icol=1,npar
	       ip=vpnt(it,icol)
	       ipar(icol)=nint(param(ip)*1e6)
               isig(icol)=nint(sigma(it)*1e6)
            enddo
	 endif
	 if (good(it)) then
	    flags=iand(int(flags),257)
	 else
	    flags=ior(int(flags),512)
	 endif

         write (20,rec=irec) jobs,satid,nxor,noalt,inc,sal,
     |                       epasc,lonasc,startt,stopt,
     |                       (ipar(i),i=1,npar),
     |                       (isig(i),i=1,npar),flags
      enddo

      write (20,rec=1) ftype,ntrack,npar,bound,iter
      close (20)

      return

 1000 i=index(xtf,' ')-1
      write (6,2000) xtf(:i)
      call fin('inimini: stop')

 2000 format ('Track catalog does not exist and I really',
     |        'absolutely need it'/'Filename : ',a)
      end
