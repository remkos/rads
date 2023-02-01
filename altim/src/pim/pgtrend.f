      subroutine pgtrend(lon,lat,c)
      include "pim.inc"
      real c(10),lat,lon

      call pgopen(device)
      call pgswin(tx(1),tx(nstack),c(9),c(10))
      call pgbox('abcnst',0.0,0,'abcnst',0.0,0)
      call pgtrend1('lat',lat ,-2.0)
      call pgtrend1('lon',lon ,-3.5)
      call pgtrend1('int',c(4),-5.0)
      call pgtrend1('slp',c(5),-6.5)
      call pgtrend1('cor',c(6),-8.0)
      call pgtrend1('fit',c(7),-9.5)
      call pgslw(3)
      call pgsci(2)
      call pgpt(nstack,tx,ty,17)
      call pgmove(tx(1),c(4))
      call pgdraw(tx(nstack),c(4)+tx(nstack)*c(5))
      call pgclos
      call pgslct(device_id)
      end

      subroutine pgtrend1(lab,x,y)
      character*80 text
      character*(*) lab
      real x,y

      write (text,600) lab,x
600   format (a,' = ',f8.3)
      call pgmtext('T',y,0.05,0.0,text)
      end
