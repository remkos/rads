c
c write Geodyn Binary I observation file to ASCII
c only for SLR data !!
c
      program slr2asc
      implicit none
      integer*2 ityp,indtim
      integer*4 isat,istat,ipre,mjd,idum1,idum2,imet,idum3
      real*4 sig,trop,cion,cmass
      real*8 fract,range
      character*80 input,output
      integer*4 iunit/10/,ounit/20/

      call getarg(1,input)
      call getarg(2,output)
      open(unit=iunit,status='old',form='unformatted',file=input)
      if (output.eq.' ' .or. output.eq.'-') then
	 ounit=6
      else
         open(unit=ounit,status='unknown',form='formatted',file=output)
      endif

  10  read (iunit,end=999) isat,ityp,indtim,istat,ipre,mjd,fract,range,
     |  idum1,idum2,sig,trop,imet,cion,idum3,cmass
      write(ounit,*) isat,ityp,indtim,istat,ipre,mjd,fract,range,
     |  idum1,idum2,sig,trop,imet,cion,idum3,cmass
      goto 10

  999 end
