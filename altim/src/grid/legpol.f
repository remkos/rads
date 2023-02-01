      subroutine legpol(ct,st,pnm,nmax)
c
c      computation of normalized Legendre polynomials for
c      colatitude teta and its derivaves to latitude:
c
c          ct = cos(teta)
c          st = sin(teta)
c
c      maximum degree equal to nmax
c      output in pnm(n,m), that is pnm(cos(teta)) = pnm(n+1,m+1)
c
c
      integer ndeg
      parameter(ndeg=360)
      integer temp,poff,coff,ioff,i,j,n,m,n2,iadr,nmax,
     |	temp1,temp2
      integer offst(0:ndeg)
      real*8 pnm(ndeg+1,ndeg+1)
      real*8 root(0:3*ndeg)
      real*8 p((ndeg+1)*(ndeg+2)/2)
      real*8 ct,st,f0,f1,f2,factor
      
      do i=0,3*ndeg
	 root(i)=dsqrt(dfloat(i))
      enddo
      
      do m=0,nmax
	 offst(m)=m*(nmax+1)-m*(m+1)/2+1
      enddo

      coff=0

      do m=0,nmax
	 ioff=offst(m)
	 if (m.eq.0) then
	    temp=ioff+m
	    p(temp)=1d0
	 else if (m.eq.1) then
	    poff=0
	    coff=ioff
	    temp=ioff+m
	    p(temp)=root(3)*st
	 else
	    poff=coff
	    coff=ioff
	    n2=m*2
	    factor=root(n2+1)/root(n2)
	    temp=poff+m-1
	    p(coff+m)=factor*st*p(temp)
	 endif
	 if (m.lt.nmax) then
	    n=m
	    iadr=ioff+n
	    factor=root(n*2+3)
	    p(iadr+1)=factor*ct*p(iadr)
	    do n=m+2,nmax
	       n2=n*2
	       f0=root(n2+1)/root(n-m)/root(n+m)
	       f1=root(n2-1)*ct
	       f2=root(n-m-1)*root(n+m-1)/root(n2-3)
	       iadr=ioff+n
	       temp1=iadr-1
	       temp2=temp1-1
	       p(iadr)=(f1*p(temp1)-f2*p(temp2))*f0
	    enddo
	 endif
      enddo

      do i=0,nmax
	 do j=0,i
	    pnm(i+1,j+1)=p(offst(j)+i)
	 enddo
      enddo
       
      return
      end
