c
c read GEODYN Trajectory File and write to ASCII format
c
      implicit none
      integer na,i,j
      character*80 input,output
      real*8 buffer(2048),buf1(48),time
      character*8 char(2000)
      logical lend

      call getarg(1,input)
      call getarg(2,output)

      open(unit=10,status='old',form='unformatted',file=input)
      open(unit=20,status='unknown',form='formatted',file=output)

c
c Read G2T header buffer
c
      read(10) buffer
      write(20,'(6d22.15)') buffer

c
c Number of alphanumeric buffers
c
      na=nint(buffer(2))

c
c Read Alphanumeric header(s)
c
      lend=.false.
      do 10 i=2,na+1
      read(10) buf1,char
      if (lend) then
         backspace(10)
         read(10) buffer
      endif
      if (lend) then
         write(20,'(6d22.15)') buffer
      else
         write(20,'(6d22.15)') buf1
         write(20,'(10a8)') char
      endif
        do 20 j=1,2000
        if (char(j)(1:6).eq.'ENDALL') lend=.true.
   20   continue
   10 continue

c
c Process Data buffers
c
      time=0
  100 read(10) buffer
      if (buffer(1).eq.9e9) goto 200
      write(20,'(6d22.15)') buffer
      goto 100

  200 continue
      write(20,'(6d22.15)') buffer
      end
