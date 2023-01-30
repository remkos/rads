c
c convert Binary Metric file to ASCII metric file
c
      implicit real*8(a-h,o-z)
      character*60 input,output
      integer*4 imout(12)
      real*8 rmout(19)

      narg=iargc()
      if (narg.ne.2) stop
     c 'usage: met2asc METRIC_BIN METRIC_ASC'

      call getarg(1,input)
      call getarg(2,output)
      open(unit=10,status='old',form='unformatted',file=input)
      open(unit=20,status='unknown',form='formatted',file=output)

c formats
  100 format(12i11)
  200 format(5d25.16)

      nobs=0
   10 read(10,end=20) imout,rmout
      write(20,100) imout
      write(20,200) rmout
      nobs=nobs+1
      goto 10

   20 write(6,*) nobs
      end
