      subroutine addhls(idx,i,j)
      integer idx,i,j,chls(3)
      real r,g,b,h,l,s
      include 'xrgb.inc'

      changed=.true.
      call grxhls(cmap(1,idx)/255.,cmap(2,idx)/255.,cmap(3,idx)/255.,
     .h,l,s)
      chls(1)=nint(h)
      chls(2)=nint(l*100)
      chls(3)=nint(s*100)
      chls(i)=chls(i)+j
      if (i.eq.1) then
         if (chls(i).ge.360) chls(i)=chls(i)-360
         if (chls(i).lt.0) chls(i)=chls(i)+360
      else
         if (chls(i).gt.100) chls(i)=100
         if (chls(i).lt.0) chls(i)=0
      endif
      call grxrgb(chls(1)*1.,chls(2)/100.,chls(3)/100.,
     .r,g,b)
      cmap(1,idx)=nint(r*255)
      cmap(2,idx)=nint(g*255)
      cmap(3,idx)=nint(b*255)
      end
