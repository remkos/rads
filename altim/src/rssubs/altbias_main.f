      integer prod,utc
      real*8  cor0,cor1,cor2,sec85,t
      logical sec
      character*80 arg
      call getarg(1,arg)
      sec=(arg(:2).eq.'-s')
10    read (5,*,end=9999) prod,t
      if (sec) then
         utc=nint(t)
      else
         utc=nint(sec85(0,utc))
      endif
      call altbias(prod,utc,cor0,cor1,cor2)
      write (*,'(i3,i14,3f12.3)')
     |		prod,utc,cor0*1d3,cor1*1d3,cor2*1d3
      goto 10
9999  end
