      implicit real*8(a-h,o-z)
      parameter(nmax=360)
      parameter(number=(nmax+1)*(nmax+2)/2)
      character*100 input,output
      character*6 char71
      real*8 cnm(number),snm(number)
      integer offset(0:nmax)

      data cnm /number*0d0/
      data snm /number*0d0/

      call getarg(1,input)
      call getarg(2,output)
      call getarg(3,out2)

      open(unit=10,status='old',form='formatted',file=input)
      open(unit=20,status='unknown',form='unformatted',file=output)

      do 10 m=0,nmax
      offset(m) = m*(nmax+1) - m*(m+1)/2 + 1
   10 continue

c read header
      read(10,*)
      read(10,*)
   20 read(10,100,end=30) char71,l,m,coefc,coefs
  100 format(a7,7x,2i3,d24.8,d15.8)
      cnm(offset(m)+l)=coefc
      snm(offset(m)+l)=coefs
      goto 20

   30 cnm(1)=1.0

c write file for 'geogen'
      write(20) (cnm(iadr),iadr=1,number)
      write(20) (snm(iadr),iadr=1,number)

      end
