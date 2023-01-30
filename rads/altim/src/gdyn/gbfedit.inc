* gbfedit.inc

      integer maxdel,ndel
      parameter (maxdel=6000)
      integer delsta(maxdel),pnt(maxdel)
      real*8 delt0(maxdel),delt1(maxdel)
      real*8 value(maxdel),val0(maxdel),val1(maxdel)
      integer*2 pnt1(10000),pnt2(10000)
      common /cdel/ delsta,pnt,delt0,delt1,value,val0,val1,
     |pnt1,pnt2,ndel
