      subroutine setcol(idx,text)
      integer idx,i,ir,ig,ib
      character*(*) text

      include 'pim.inc'

* Set colour index idx based on ASCII RGB colour representation
* Examples of input (text):
* - RGB values (0-255)
*   255 255 255
* - RGB values (Hexadecimal):
*   #ffffff
* - HLS values (0-360, 0-100, 0-100):
*   180 100 100 HLS

      i=index(text,'#')
      if (i.gt.0) then
         read (text(i+1:i+6),'(3z2)') ir,ig,ib
      else
         read (text,*,iostat=ios) ir,ig,ib
      endif
      if (index(text,'HLS').gt.0) then
         call grxrgb(ir*1.,ig/100.,ib/100.,rgb(1,i),rgb(2,i),rgb(3,i))
      else
         rgb(1,idx)=ir/255.
         rgb(2,idx)=ig/255.
         rgb(3,idx)=ib/255.
      endif
      call pgscr(idx,rgb(1,idx),rgb(2,idx),rgb(3,idx))
      end
