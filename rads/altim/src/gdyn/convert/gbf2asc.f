c
c     read GEODYN II binary format and convert to ASCII
c
      implicit real*8(a-h,o-z)
      real*8 buf(200,10)
      character*60 input,output

      equivalence (rprepro,i4)

      call getarg(1,input)
      call getarg(2,output)

      open(unit=10,status='old',form='unformatted',file=input)
      open(unit=20,status='unknown',form='formatted',file=output)

c     read and write buffers
   10 read(10,end=999) buf
      write(20,200) buf
  200 format(6d24.17)
      goto 10

  999 end
