      subroutine annotate(rlon,rlat,z,mx,my,dc,ic)
      include "pim.inc"
      integer mx,my,ic
      real*4  z(mx,my),dx,dy,rlon,rlat
      integer ipower/3/,kx,ky,n
      character text*10
      real pi,rad
      real xg,yg,rx,ry
      real z00,z01,z10,z11,dzdx,dzdy,angle,cosa,sina
      real zi,dc,bx,by,xbox(4),ybox(4)
      real xp0,xp1,yp0,yp1
      real xblc,xtrc,yblc,ytrc,xsp,ysp
      
      call pgqvp(3,xp0,xp1,yp0,yp1)
      call pgqwin(xblc,xtrc,yblc,ytrc)
      call grchsz(1,xsp,ysp,xsp,ysp)
      pi=4*atan(1.)
      rad=pi/180
      dx=(xw1-xw0)/(mx-1)
      dy=(yw1-yw0)/(my-1)
      xg=(rlon-xw0)/dx+1.
      yg=(rlat-yw0)/dy+1.
      kx=xg
      ky=yg
      rx=xg-kx
      ry=yg-ky
      z00=z(kx,ky)
      z01=z(kx,ky+1)
      z10=z(kx+1,ky)
      z11=z(kx+1,ky+1)
      dzdx=z10-z00+z11-z01
      dzdy=z01-z00+z11-z10
      angle=atan2(dzdy,-dzdx)/rad
      call pmcvec(1,rlon,rlat,angle)
      if (angle.gt.90 .and. angle.lt.270) angle=angle-180
      if (angle.gt.-270 .and. angle.lt.-90) angle=angle+180
      zi=(1-rx)*((1-ry)*z00+ry*z01)+rx*((1-ry)*z10+ry*z11)
      zi=nint((zi-rminc)/dc)*dc+rminc
      call pgnumb(nint(zi*(10**ipower)),-ipower,1,text,n)
      bx=.5*xsp/(xp1-xp0)*(xtrc-xblc)*n
      by=.5*ysp/(yp1-yp0)*(ytrc-yblc)
      sina=sin(angle*rad)
      cosa=cos(angle*rad)
      xbox(1)=rlon-bx*cosa+by*sina
      ybox(1)=rlat-bx*sina-by*cosa
      xbox(2)=rlon+bx*cosa+by*sina
      ybox(2)=rlat+bx*sina-by*cosa
      xbox(3)=rlon+bx*cosa-by*sina
      ybox(3)=rlat+bx*sina+by*cosa
      xbox(4)=rlon-bx*cosa-by*sina
      ybox(4)=rlat-bx*sina+by*cosa
      call pgsci(c_bg)
      call pgpoly(4,xbox,ybox)
      call pgsci(ic)
      call pgptxt(rlon+.8*by*sina,rlat-.8*by*cosa,angle,.5,text(:n))
      end
