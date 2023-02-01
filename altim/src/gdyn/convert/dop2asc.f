c
c write Geodyn Binary I observation file to ASCII
c only for DOPPLER data !!
c
      implicit real*8 (a-h,o-z)
      real*4 sig,trop,cion,cmass
      integer*2 ityp,indtim
      character*60 input,output

      call getarg(1,input)
      call getarg(2,output)

      open ( 10 , file = input , form = 'unformatted' )
      open ( 20 , file = output , form = 'formatted' )

      nrec=0
  10  read (10,end=100) isat,ityp,indtim,istat,ipre,mjd,fract,rrate,
     .  idum1,idum2,sig,icount,trop,cion,idum3,cmass
      write(20,*) isat,ityp,indtim,istat,ipre,mjd,fract,rrate,
     .  idum1,idum2,sig,icount,trop,cion,idum3,cmass
      nrec=nrec+1
      goto 10

  100 write(6,*) nrec
      end
