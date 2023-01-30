c
c     read GEODYN II observation partail derivative file
c     and convert from Binary to ASCII
c
      implicit real*8(a-h,o-z)
      parameter(nmxp=10000)
      real*8 header(30),rgroup(100),rnrgrp(100)
      real*8 apriori(nmxp),value(nmxp),label(3*nmxp)
      real*8 rinfo(50),partial(nmxp)
      real*8 rms(100),rmn(100)
      integer nrobs(100)
      character*60 input,out
      character*120 text

      data rms /100*0d0/
      data rmn /100*0d0/
      data nrobs /100*0/

      nargc=iargc()
      if (nargc.ne.2) stop
     c  'usage: parfil2asc PARTIALS.GEODYN PARTIALS.ASC'

      call getarg(1,input)
      call getarg(2,out)

      open(unit=10,status='old',form='unformatted',file=input)
      open(unit=20,status='unknown',form='formatted',file=out)

c Formats
  200 format(6d24.17)

c Header record
      read(10) header
      write(20,200) header
c File number
      rfilnr=header(1)
c Ematrx number
      efilnr=header(2)
c Number of estimated parameters
      nparam=nint(header(3))
c Number of weighted observations
      nobs=nint(header(5))
c Total variance
      var=header(6)
c Weighted variance
      wvar=header(7)
c Arc variance
      arcvar=header(8)
c Number of satellites
      nsat=nint(header(9))
c Central forcing body
      ibody=nint(header(10))
c Date of job
      datej=header(11)
c Time of job
      timej=header(12)
c Geodyn version
      iversion=nint(header(14))
c Number of parameter groups
      ngi=nint(header(15))
c Speed of light
      vlight=header(16)

c Parameter group identifiers
      read(10) (rgroup(i),i=1,ngi)
      write(20,200) (rgroup(i),i=1,ngi)

c Number of parameter records
      read(10) (rnrgrp(i),i=1,ngi)
      write(20,200) (rnrgrp(i),i=1,ngi)

c Parameter labels record
      read(10) rnr
      backspace(10)
      nr=nint(rnr/3)
      read(10) rnr,(label(i),i=1,3*nr)
      write(20,200) rnr,(label(i),i=1,3*nr)

c Parameter values record
      read(10) rnr
      backspace(10)
      nr=nint(rnr)
      read(10) rnr,(value(i),i=1,nr)
      write(20,200) rnr,(value(i),i=1,nr)

c Parameter values record
      read(10) rnr
      backspace(10)
      nr=nint(rnr)
      read(10) rnr,(apriori(i),i=1,nr)
      write(20,200) rnr,(apriori(i),i=1,nr)

c Measurement partial derivative records
      nrec=0
   60 read(10,end=999) rinfo,(partial(i),i=1,nparam)
      write(20,200) rinfo,(partial(i),i=1,nparam)
      if (dabs(rinfo(1)+999999.0).lt.0.1) then
         write(6,*) 'NORMAL TERMINATION: SENTINEL RECORD REACHED'
         write(6,*) 'NUMBER OF RECORDS: ',nrec
         write(6,*)
         goto 999
      endif
      nrec=nrec+1
      if (mod(nrec,1000).eq.0) write(6,*) nrec
c observation time
      tobs=rinfo(1)+rinfo(2)
c observation type
      mtype=int(rinfo(5))
c observation value
      val=rinfo(13)
c sum of observation time corrections
      timcor=rinfo(14)
c sum of observation corrections
      obscor=rinfo(15)
c observation residual
      resid=rinfo(17)
c observation sigma
      if (rinfo(18).eq.0) then
         write(6,*) 'This measurement not used'
         write(6,*) 'Measurement type and residual: ',mtype,resid
         write(6,*)
      else
         sigma=dsqrt(1d0/rinfo(18))
         nrobs(mtype)=nrobs(mtype)+1
         rmn(mtype)=rmn(mtype)+resid
         rms(mtype)=rms(mtype)+resid**2
      endif
      goto 60

  999 write(6,*) 'Residual statistics'
      write(6,*)
      do 70 i=1,100
      if (mod(i,2).eq.1) fac=1d2
      if (mod(i,2).eq.0) fac=1d3
      if (nrobs(i).ne.0) then
         write(6,'(i4,i6,2f15.4)') i,nrobs(i),fac*rmn(i)/nrobs(i),
     c                fac*dsqrt(rms(i)/nrobs(i))
      endif
   70 continue

      end
