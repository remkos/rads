c
c
c    E02ADF
c
c    Purpose: approximate a set of (x,y) points by a polynomial
c             of wanted degree using the least squares method
c
c    Usage:   call e02adf(x,y,m,k1,sig,p,c)
c
c    Parameters:
c        x  : i x(m) contains the abscissa of the given points
c        y  : i y(m) contains the ordinates of the given points
c        m  : i the number of points
c        k1 : i the wanted degree of the polynomial + 1
c        sig: o sig(k1), see under remarks
c        p  : o p(k1) contains the computed coefficients
c               of the polynomial sum (j:1,k1) p(j)*x**(j-1)
c        c  :   working area of dimension of at least 6+6*k+2*m
c
c    Remarks:
c        array sig(k1) contains the following quantities:
c         sig(1)=d1/(m-1), sig(i)=di/(m-i-1) i=2(1)k1,
c        where di= sum(j:1,m) (y(j)-fi(x(j)))**2
c        and where fi(x) is a polynomial of degree i-1
c
c    References:
c        RC-TWA-75003 NUMLIBDA Methodebeschrijvingen
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine e02adf(x,y,m,k1,sig,p,c)
      implicit double precision (a-h,o-z)
      logical swx
      dimension x(m),y(m),sig(k1),p(k1),c(*)
      data z,u/0.d0,1.d0/
      k=k1-1
      ic=k1+2
      iz=ic+k1
      ib=iz+k1
      is=ib+k1
      ia=is+k
      i1=ia+k
      i2=i1+m
      swx=.true.
      c(ib)=z
      c(1)=z
      c(2)=z
      del=z
      ome=z
      c(ic)=u
      r2=m
      do 10 i=1,m
      del=del+y(i)*y(i)
      c(i+i1)=z
      c(i+i2)=u
   10 ome=ome+y(i)
      c(is)=ome/r2
      c(iz)=c(is)
      del=del-c(is)*ome
      sig(1)=del/(m-u)
      do 80 i=0,k-1
      d=z
      do 20 j=1,m
   20 d=d+x(j)*c(j+i2)**2
      c(i+1+ia)=d/r2
      r1=r2
      r2=z
      ome=z
      ru=z
      do 30 j=1,m
      d=c(i+ib)*c(j+i1)
      c(j+i1)=c(j+i2)
      c(j+i2)=(x(j)-c(i+1+ia))*c(j+i2)-d
      r2=r2+c(j+i2)**2
      ome=ome+y(j)*c(j+i2)
   30 ru=ru+x(j)*c(j+i2)*c(j+i1)
      c(i+1+ib)=ru/r1
      c(i+1+is)=ome/r2
      del=del-c(i+1+is)*ome
      sig(i+2)=del/(m-i-1)
      if ((sig(i+2).lt.sig(i+1)).and.(swx)) goto 60
      if (swx) goto 40
      if (sig(i+2).le.g) swx=.true.
      goto 60
   40 continue
      do 50 j=0,i
   50 p(j+1)=c(j+iz)
      swx=.false.
      g=sig(i+1)
      n=i
   60 continue
      do 70 j=0,i
      d=c(j+2)*c(i+ib)
      c(j+2)=c(j+ic)
      c(j+ic)=c(j+1)-c(i+1+ia)*c(j+ic)-d
   70 c(j+iz)=c(j+iz)+c(i+1+is)*c(j+ic)
      c(i+1+iz)=c(i+1+is)
      c(i+1+ic)=u
      c(i+3)=z
   80 continue
      if (swx) then
      do 90 j=0,k
   90 p(j+1)=c(j+iz)
      n=k
      endif
      do 100 i=n+1,k
  100 p(i+1)=z
      return
      end
c
c
c    E02BAF
c
c    Purpose
c       to interpolate, differentiate and integrate a tabulated
c       function on non-equidistand points by cubic splines
c
c    Usage
c       call e02baf(n,m,g,a,x,z,f,f1,v1,vn,c,d)
c
c    Description of parameters
c       n    i m=-1 given points are x1....xn, n>5
c            i m>=1 upper boundery of quadrature is xn
c       m    i m=-1 array a is filled, must occur once
c            i m<-1 f(x) and f'(x) are calculated
c            i m>=1 integration on interval (xm,xn)
c       g    i array g(n) contains the functionvalues
c       a    o array a(n) contains the moments at points xk
c              and is filled if (m.eq.-1), and is used when (m.ne.-1)
c       x    i x(n) contains the abcissae for which f(x) and f'(x)
c              are calculated if m<-1
c       z    i if m<-1 f(z) and f'(z) are calculated
c       f    o m<-1 value f(z)
c            o m>=1 value of quadrature
c       f1   o m<-1 value f'(z)
c       v1     v1(n) auxiliary workingspace
c       vn     vn(n) auxiliary workingspace
c       c      c(n)  auxiliary workingspace
c       d      d(n)  auxiliary workingspace
c
c    E02BAF uses subroutine SPICU1
c
c    Remarks
c       if spicun is used, then first the array a must be filled
c       with the moments. this is done by executing spicun with
c       m=-1 once. after this spicun can be used for interpolation
c       differentiation and intergration.
c
c    Method
c       refer to rc-twa-76001 of computer centre d.u.t.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine e02baf(n,m,g,a,x,z,f,f1,v1,vn,c,d)
      implicit double precision (a-h,o-z)
      dimension g(n),a(n),x(n),v1(n),vn(n),c(n),d(n)
      data zz,uu/0.d0,1.d0/
      n1=n-1
      if (m.eq.-1) then
      do 10 i=2,n
   10 a(i)=6.d0*(g(i)-g(i-1))/(x(i)-x(i-1))
      do 20 i=2,n1
   20 a(i)=a(i+1)-a(i)
      t=2.d0*(x(3)-x(1))
      d(2)=uu/t
      do 30 i=2,n-2
      c(i)=(x(i+1)-x(i))*d(i)
      t=2.d0*(x(i+2)-x(i))-t*c(i)*c(i)
   30 d(i+1)=uu/t
      do 40 i=2,n1
      v1(i)=zz
   40 vn(i)=zz
      v1(2)=uu
      vn(n1)=uu
      call spicu1(n,1,v1,c,d)
      call spicu1(n,0,vn,c,d)
      call spicu1(n,1,a,c,d)
      r=zz
      t=zz
      s=zz
      v=zz
      p=zz
      q=zz
      do 70 i=2,5
      u=uu
      w=uu
      do 50 j=2,5
      if (j.ne.i) then
      u=u*(x(1)-x(j))
      w=w*(x(i)-x(j))
      endif
   50 continue
      u=u/w
      r=r+u*v1(i)
      v=v+u*vn(i)
      p=p+u*a(i)
      u=uu
      w=uu
      do 60 j=1,4
      if (j.ne.i-1) then
      u=u*(x(n)-x(n-j))
      w=w*(x(n-i+1)-x(n-j))
      endif
   60 continue
      u=u/w
      s=s+u*v1(n-i+1)
      t=t+u*vn(n-i+1)
      q=q+u*a(n-i+1)
   70 continue
      w=x(2)-x(1)
      u=x(n)-x(n1)
      r=uu+w*r
      v=u*v
      s=w*s
      t=uu+u*t
      e=r*t-s*v
      a(1)=(t*p-q*v)/e
      a(n)=(r*q-p*s)/e
      r=w*a(1)
      s=u*a(n)
      do 80 i=2,n1
   80 a(i)=a(i)-r*v1(i)-s*vn(i)
      else
      if (m.lt.-1) then
      do 90 i=n,1,-1
      k=i
      if (z.ge.x(k)) goto 100
   90 continue
  100 i=k+1
      w=x(i)-x(k)
      t=(z-x(k))/w
      s=t*(uu-t)*w*w/6.d0
      r=a(i)*s
      s=a(k)*s
      f=t*(g(i)-r)+(uu-t)*(g(k)-s)-(r+s)
      s=t*t*a(i)-(uu-t)*(uu-t)*a(k)-(a(i)-a(k))/3.d0
      f1=(g(i)-g(k))/w+0.5d0*w*s
      else
      s=zz
      t=zz
      do 110 i=m+1,n
      k=i-1
      w=x(i)-x(k)
      s=s+(a(i)+a(k))*w*w*w
  110 t=t+(g(i)+g(k))*w
      f=0.5d0*t-s/24.d0
      endif
      endif
      return
      end
      subroutine spicu1(n,l,b,c,d)
      double precision b,c,d
      dimension b(n),c(n),d(n)
      n1=n-1
      if (l.eq.1) then
      do 10 i=3,n1
   10 b(i)=b(i)-c(i-1)*b(i-1)
      do 20 i=2,n1
   20 b(i)=b(i)*d(i)
      else
      b(n1)=d(n1)
      endif
      do 30 i=n-2,2,-1
   30 b(i)=b(i)-c(i)*b(i+1)
      return
      end
