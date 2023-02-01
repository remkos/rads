      program lgbf
      integer*2 j5(2)
      integer*4 idata(17),j1,j9,j13,j17,j37,j41,j53
      real*4    std,trop,iono,displ,cmas
      real*8	tfrac,range
      equivalence
     |  (idata(1),j1),    (idata(2),j5(1)),  (idata(3),j9),
     |  (idata(4),j13),   (idata(5),j17),    (idata(6),tfrac),
     |  (idata(8),range), (idata(10),j37),   (idata(11),j41),
     |  (idata(12),std),  (idata(13),trop),  (idata(14),j53),
     |  (idata(15),iono), (idata(16),displ), (idata(17),cmas)
      character*80 arg
      
      call getarg(1,arg)
      open (10,file=arg,status='old',form='unformatted')
10    read (10,end=9999) idata
      write (*,*) j1,j5,j9,j13,j17,tfrac,range,j37,j41,std,trop,
     |j53,iono,displ,cmas
      goto 10
9999  end
