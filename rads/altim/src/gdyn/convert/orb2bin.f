c
c read ASCII GEODYN Trajectory File and write to Binary
c
      implicit none
      integer na,i,j
      character*80 input,output
      real*8 buffer(2048),buf1(48),time
      character*8 char(2000)
      logical lend

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
      lend=.false.
      do 10 i=2,na+1
      if (lend) then
         read(10,'(6d22.15)') buffer
         write(20) buffer
      else
         read(10,'(6d22.15)') buf1
         read(10,'(10a8)') char
         write(20) buf1,char
      endif
         do 20 j=1,2000
         if (char(j)(1:6).eq.'ENDALL') lend=.true.
   20    continue
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
