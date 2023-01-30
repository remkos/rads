      program xtfrewr
      implicit none
      include "data.inc"
      include "grid.inc"
      include "stat.inc"
      integer i,irec

      character	ftype*4,new*80
      integer*4 inc,sal,epasc,lonasc,startt,stopt,
     |          ipar(maxpar),isig(maxpar),mpar
      integer*2 jobs,satid,nxor,noalt,flags
      
      call getarg (1,xtf)
      read (xtf,*) mpar
      call getarg (2,xtf)
      call getarg (3,new)

* Open track catalog

      open (20,file=xtf,form='unformatted',access='direct',
     |      recl=58,status='old')
      read (20,rec=1) ftype,ntrack,npar

      if (ftype.eq.'@XTB') then
         short=.true.
	 xaflen=28
         xxflen=48
         xtflen=34+npar*8
         close (20)
         open (20,file=xtf,form='unformatted',access='direct',
     |         recl=xtflen,status='old')
         read (20,rec=1) ftype,ntrack,npar,bound,iter
         open (30,file=new,form='unformatted',access='direct',
     |      recl=34+mpar*8,status='new')
      else
	 call fin('inimini: only @XTB supported')
      endif

* Start reading catalog

      do irec=2,ntrack+1
         read (20,rec=irec) jobs,satid,nxor,noalt,inc,sal,
     |                      epasc,lonasc,startt,stopt,
     |                      (ipar(i),i=1,npar),
     |                      (isig(i),i=1,npar),flags
         write (30,rec=irec) jobs,satid,nxor,noalt,inc,sal,
     |                       epasc,lonasc,startt,stopt,
     |                       (ipar(i),i=1,mpar),
     |                       (isig(i),i=1,mpar),flags
      enddo

      write (30,rec=1) ftype,ntrack,mpar,bound,iter
      close (20)
      close (30)

      end
