c
c     Conversion of binary "matrix" file to ASCII
c
      implicit real*8(a-h,o-z)
      parameter(nparam=20000,nmxsat=50)
      character*60 input,output

      real*8 header(30),label(3,nparam)
      real*8 tiieout(nmxsat*(6*nparam+6))
      real*8 row(nparam)

      call getarg(1,input)
      call getarg(2,output)

      open(unit=10,status='old',form='unformatted',file=input)
      open(unit=20,status='unknown',form='formatted',file=output)

c     header record
      read(10) header
      number=nint(header(3))
      ntiieout=nint(header(19))
      write(6,*) number
      if (number.gt.nparam) stop
     c 'increase number size'
      write(20,'(6d22.15 )') (header(i),i=1,2)
      write(20,'(6d22.15)') (header(i),i=3,30)

c     parameter group identifier
      read(10) (row(i),i=1,nint(header(15)))
      write(20,'(6d22.15)') (row(i),i=1,nint(header(15)))

c     number of parameters record
      read(10) (row(i),i=1,nint(header(15)))
      write(20,'(6d22.15)') (row(i),i=1,nint(header(15)))

c     labels
      if (ntiieout.gt.0) then
         read(10) dword,(label(1,j),j=1,number),
     c    (label(2,j),j=1,number),(label(3,j),j=1,number),dresid,
     c    tiieout1,tiieout2,rnsat
          nsat=nint(rnsat)
          if (nsat.gt.nmxsat) stop
     c      'increase maximum number of satellites'
          nhelp=nsat*(6*number+6)
          backspace(10)
          read(10) dword,(label(1,j),j=1,number),
     c    (label(2,j),j=1,number),(label(3,j),j=1,number),dresid,
     c    tiieout1,tiieout2,rnsat,rgm,ridsat,(tiieout(i),i=1,nhelp)
c for preventing underflows
          do 15 i=1,nhelp
             if (dabs(tiieout(i)).lt.1d-100) tiieout(i)=0d0
   15     continue
          write(20,'(6d22.15)') dword,(label(1,j),j=1,number),
     c    (label(2,j),j=1,number),(label(3,j),j=1,number),dresid,
     c    tiieout1,tiieout2,rnsat,rgm,ridsat,(tiieout(i),i=1,nhelp)
      else
         read(10) dword,(label(1,j),j=1,number),
     c    (label(2,j),j=1,number),(label(3,j),j=1,number),dresid
         write(20,'(6d22.15)') dword,(label(1,j),j=1,number),
     c   (label(2,j),j=1,number),(label(3,j),j=1,number),dresid
      endif

c     parameters values record (=a priori)
      read(10) dnum,(row(i),i=1,number)
      write(20,'(6d22.15)') dnum,(row(i),i=1,number)

c     parameter a priori values record
      read(10) dnum,(row(i),i=1,number)
      write(20,'(6d22.15)') dnum,(row(i),i=1,number)

c     normal equations
      do 10 i=1,number
      read(10) drow,delem,(row(j),j=i,number),resid
      write(20,'(6d22.15)') drow,delem,(row(j),j=i,number),resid
      if (mod(i,100).eq.0) write(6,*) i
   10 continue

      end

