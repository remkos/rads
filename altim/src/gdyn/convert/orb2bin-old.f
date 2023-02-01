c
c read ASCII GEODYN Trajectory File and write to Binary
c
      implicit real*8(a-h,o-z)
      character*60 input,output
      real*8 buffer(2048)
      character*8 char(2046)

      call getarg(1,input)
      call getarg(2,output)

      open(unit=10,status='old',form='formatted',file=input)
      open(unit=20,status='unknown',form='unformatted',file=output)

c
c Read G2T header buffer
c
      read(10,'(6d22.15)') buffer
      write(20) buffer

c
c Number of alphanumeric buffers
c
      na=nint(buffer(2))

c
c Read Alphanumeric header(s)
c
      do 10 i=2,na+1
      read(10,'(2d22.15)') dum1,dum2
      read(10,'(10a8)') char
      write(20) dum1,dum2,char
   10 continue

c
c Process Data buffers
c
      time=0
  100 read(10,'(6d22.15)') buffer
      if (buffer(1).eq.9d9) goto 200
      write(20) buffer
      goto 100

  200 continue
      write(20) buffer
      end
