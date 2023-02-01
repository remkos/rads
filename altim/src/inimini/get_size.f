      subroutine get_size
      implicit none
      include "stat.inc"
      include "grid.inc"
      include "data.inc"
      integer l
      character*4 chriter,ftype,error*132
      equivalence (iter,chriter)

* First open catalog

      l=index(xtf,' ')-1
      open (20,file=xtf,form='unformatted',access='direct',
     |      recl=58,status='old',err=1010)
      read (20,rec=1) ftype,ntrack,npar

      if (ntrack.gt.mtrack)call fin('inimini: too many tracks')

* Check filetype, file must be reopened with proper length
* if it is a variable length file (extended file)

      if (ftype.eq.'@XTB') then
         short=.true.
	 xaflen=28
         xxflen=48
         xtflen=34+npar*8
         close (20)
         open (20,file=xtf,form='unformatted',access='direct',
     |         recl=xtflen,status='old')
         read (20,rec=1,err=1020) ftype,ntrack,npar,bound,iter
      else if (ftype.eq.'@XTF') then
	 call fin('inimini: @XTF no longer supported')
         extend=.false.
         open (20,file=xtf,form='unformatted',access='direct',
     |         recl=58,status='old')
         read (20,rec=1,err=1020) ftype,ntrack
         npar=2
      else if (ftype.eq.'@XTE') then
         extend=.true.
	 xaflen=28+(npar-1)*4
         xtflen=34+npar*8
         xxflen=42+(npar-1)*8
         close (20)
         open (20,file=xtf,form='unformatted',access='direct',
     |         recl=xtflen,status='old')
         read (20,rec=1,err=1020) ftype,ntrack,npar,bound,iter
      else
         goto 1000
      endif

      close(20)

      if (chriter.eq.'    ' .or. init) iter=0

      if (minlon.eq.maxlon) then
         minlon=bound(1)
         maxlon=bound(2)
      else
	 bound(1)=nint(minlon)
	 bound(2)=nint(maxlon)
      endif
      if (minlat.eq.maxlat) then
         minlat=bound(3)
         maxlat=bound(4)
      else
	 bound(3)=nint(minlat)
	 bound(4)=nint(maxlat)
      endif

      return

 1000 write (error,1100) xtf(:l)
      call fin(error)
 1010 write (error,1110) xtf(:l)
      call fin(error)
 1020 write (error,1120) xtf(:l)
      call fin(error)
 1100 format ('inimini: file ',a,' is not in XTB or XTE format')
 1110 format ('inimini: error opening XTF/XTE file: ',a)
 1120 format ('inimini: error reading XTF/XTE file: ',a)
      end
