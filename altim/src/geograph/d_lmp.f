*+D_LMP -- Initialize constant matrix Dlmp
**
      SUBROUTINE D_LMP
*
* Compute F_lmp and then build up the matrix Dlmp, as defined in
* Rosborough (5.3).
* The definition by Wagner (1985) for H_lmk is applied.
* Also some other workspaces are filled.
*-
* 18-Mar-1992 - Created. Remko Scharroo.
*  3-Aug-2001 - More comments.
*------------------------------------------------------------------------
      include "geograph.inc"
      real*8 alpha,beta,fac,s
      integer l,m,p,i,k,m1

* Generate the lookup table for F_lmp. Store it temporarily in array dlmp.

      if (test) write (*,*) 'Initialize Flmp and Dlmp'
      call f_lmp(lmax,incl,pnt,dlmp)

* Now first fill some workspaces:
* cosi(k) = cos(incl)**k  for  k=0..lmax
* sini(k) = sin(incl)**k  for  k=0..lmax*2
* bin(m,k) = binomial coefficient (m,k)  for  m=0..lmax,k=0..m
* h(l) = l*(l-1)/2  for  l=0..lmax

      cosi(0)=1
      s=cos(incl)
      do k=1,lmax
         cosi(k)=cosi(k-1)*s
      enddo

      sini(0)=1
      s=sin(incl)
      do k=1,2*lmax
         sini(k)=sini(k-1)*s
      enddo

      do m=0,lmax
         bin(m,0)=1
         bin(m,m)=1
         m1=m+1
         do k=1,m/2
            bin(m,k)=bin(m,k-1)*(m1-k)/k
            bin(m,m-k)=bin(m,k)
         enddo
      enddo

      h(0)=1
      do l=1,lmax
         h(l)=h(l-1)+l
      enddo

* Build up Dlmp.
* When the perturbation frequency is close to -1, 0, or +1 cpr, Dlmp is
* set to 0. What is close is determined by the parameter deadband in
* cpr.

      do l=1,lmax
         fac=a0*(ae/a0)**l
         if (test) write (*,620) l,a0,ae,n0,ogdot,fac
         do m=0,l
            do p=0,l
               k=l-2*p
               beta=(k*wmdot+m*ogdot)/n0
               if (dabs(beta-1).le.deadband) then
                  alpha=0d0
               else if (dabs(beta+1).le.deadband) then
                  alpha=0d0
               else if (dabs(beta).le.deadband) then
                  alpha=0d0
               else
                  alpha=(l+1-2*k/beta)/(beta*beta-1)
               endif
               i=pnt(l)+m*(l+1)+p
               if (test) write (*,600) l,m,p,k,
     |            beta,alpha,dlmp(i),fac*dlmp(i)*alpha
               dlmp(i)=fac*dlmp(i)*alpha
            enddo
         enddo
      enddo

  600 format (3i3,i4,2f15.9,f15.9,f15.3)
  620 format ('l,a0,ae,n0,ogdot,fac  :',i3,2f14.3,2f15.12,f20.9/
     |'  l  m  p   k',10x,' beta',10x,'alpha',10x,'F_lmp',10x,'D_lmp')
      end
