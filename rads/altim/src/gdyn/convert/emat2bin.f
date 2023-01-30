c
c     Conversion of ASCII "matrix" file to binary
c
      implicit real*8(a-h,o-z)
      parameter(nparam=20000)
      character*60 input,output

      real*8 header(30),label(3,nparam)
      real*8 row(nparam)

      call getarg(1,input)
      call getarg(2,output)

      open(unit=10,status='old',form='formatted',file=input)
      open(unit=20,status='new',form='unformatted',file=output)

c     header record
      read(10,'(6d15.8 )') (header(i),i=1,2)
      read(10,'(6d22.15)') (header(i),i=3,30)
      number=nint(header(3))
      write(6,*) number
      if (number.gt.nparam) stop
     c 'increase number size'
      write(20) header

c     parameter group identifier
      read(10,'(6d22.15)') (row(i),i=1,nint(header(15)))
      write(20) (row(i),i=1,nint(header(15)))

c     number of parameters record
      read(10,'(6d22.15)') (row(i),i=1,nint(header(15)))
      write(20) (row(i),i=1,nint(header(15)))

c     labels
      read(10,'(6d22.15)') dword,(label(1,j),j=1,number),
     c (label(2,j),j=1,number),(label(3,j),j=1,number),dresid
      write(20) dword,(label(1,j),j=1,number),
     c (label(2,j),j=1,number),(label(3,j),j=1,number),dresid

c     parameters values record (=a priori)
      read(10,'(6d22.15)') dnum,(row(i),i=1,number)
      write(20) dnum,(row(i),i=1,number)

c     parameter a priori values record
      read(10,'(6d22.15)') dnum,(row(i),i=1,number)
      write(20) dnum,(row(i),i=1,number)

c     normal equations
      do 10 i=1,number
      read(10,'(6d22.15)') drow,delem,(row(j),j=i,number),resid
      write(20) drow,delem,(row(j),j=i,number),resid
      if (mod(i,100).eq.0) write(6,*) i
   10 continue

      end

