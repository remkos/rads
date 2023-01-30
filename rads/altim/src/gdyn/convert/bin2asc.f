      implicit real*8(a-h,o-z)
      real*8 x(6),gtrack(6)
      character*100 input,output

      call getarg(1,input)
      call getarg(2,output)

      open(unit=10,status='old',form='unformatted',file=input)
      open(unit=20,status='unknown',form='formatted',file=output)

   10 read (10,end=20) idate,rsec,x,gtrack
      write(20,200) idate,rsec,x,gtrack
  200 format(i6,f10.1,12d17.10)
      goto 10

   20 end
