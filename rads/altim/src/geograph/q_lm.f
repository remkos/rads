*+Q_LM - Initialize matrix Qlm
**
      SUBROUTINE Q_LM (PHI)
      REAL*8 PHI
*
* Initialize the latitude dependent matrix Q_lm, as given by Rosborough
* (5.44-45). This requires also the computation of Y_c, Y_s, Psi_c,
* Psi_s, Phi_c and Phi_s as defined by (5.36-41).
*-
* 19-Mar-1993 - Created by Remko Scharroo
*  3-Aug-2001 - Added more comments
*  9-Aug-2001 - Combined q_lm with y_psi and phi_lmp
*----------------------------------------------------------------------
      include "geograph.inc"
      real*8  c,s,u,v,wk1(0:ndeg),wk2(0:ndeg),
     |        phic(-ndeg:ndeg),phis(-ndeg:ndeg),
     |        yc(0:ndeg),ys(0:ndeg),psic(0:ndeg),psis(0:ndeg)
      integer i,k,l,m

* Compute Y_c, Y_s, Psi_c and Psi_s but first fill workspaces:
*     wk1(k) = (sin(i)**2-sin(phi)**2)**(k/2)
*   wk2(2*k) = (-1)**k*sin(phi)**(2*k)
* wk2(2*k+1) = (-1)**k*sin(phi)**(2*k+1)

      s=sin(phi)
      c=sqrt(sini(2)-s**2)
      wk1(0)=1d0
      do k=1,lmax
         wk1(k)=wk1(k-1)*c
      enddo

      c=-s**2
      wk2(0)=1
      do k=2,lmax,2
         wk2(k)=wk2(k-2)*c
      enddo

      wk2(1)=s
      do k=3,lmax,2
         wk2(k)=wk2(k-2)*c
      enddo

      do m=0,lmax

* Compute Y_c and Psi_c (i=2*k)

         u=0d0
         v=0d0
         do i=0,m,2
            s=bin(m,i)*wk1(m-i)*wk2(i)
            u=u+s
            v=v+s*cosi(i)
         enddo
         yc(m)=u
         psic(m)=v

* Compute Y_s and Psi_s (i=2*k+1)

         u=0d0
         v=0d0
         do i=1,m,2
            s=bin(m,i)*wk1(m-i)*wk2(i)
            u=u+s
            v=v+s*cosi(i)
         enddo
         ys(m)=u
         psis(m)=v
      enddo

* Print result for Y_c, Y_s, Psi_c and Psi_s when requested

      if (test) then
         write (*,600)
         do m=0,lmax
            write (*,610) m,yc(m),ys(m),psic(m),psis(m)
         enddo
      endif

* Fill workspace: wk1(k) = cos(phi)**k.

      wk1(0)=1d0
      s=cos(phi)
      do k=1,lmax
         wk1(k)=wk1(k-1)*s
      enddo
      
      do m=0,lmax

* Compute Phi_c and Phi_s, as defined by (5.36-37) of Rosborough.
* In the formulations of Phi_c and Phi_s is the term (l-2p)/|l-2p|.
* Here K takes the place of kappa=|l-2p|.
* We build up Phi_c and Phi_s for a single value of M and store them for
* positive and negative K as -K and +K.

         do k=0,lmax
            s=wk1(m)*sini(m+k)
            phic(+k)=(+yc(k)*psic(m)+ys(k)*psis(m))/s
            phic(-k)=(+yc(k)*psic(m)-ys(k)*psis(m))/s
            phis(+k)=(+ys(k)*psic(m)-yc(k)*psis(m))/s
            phis(-k)=(-ys(k)*psic(m)-yc(k)*psis(m))/s
         enddo

* Print results for Phi_c and Phi_s when requested

         if (test) then
            write (*,620)
            do k=0,lmax
               write (*,630) m,k,phic(k),phic(-k),phis(k),phis(-k)
            enddo
         endif

* Finally compute Q_c and Q_s
      
         do l=max(1,m),lmax
            i=pnt(l)+m*(l+1)
            c=0d0
            s=0d0
            do k=l,-l,-2
               c=c+dlmp(i)*phic(k)
               s=s+dlmp(i)*phis(k)
               i=i+1
            enddo
            i=h(l)+m
            if (mod(l-m,2).eq.0) then
               qc(i)=c
               qs(i)=-s
            else
               qc(i)=s
               qs(i)=c
            endif
         enddo

* Print results for Q_c and Q_s when requested

         if (test) then
            write (*,640)
            do l=max(1,m),lmax
               i=h(l)+m
               write (*,650) l,m,qc(i),qs(i)
            enddo
         endif
      enddo

600   format ('  m',9x,'y^c(m)',9x,'y^s(m)',7x,'Psi^c(m)',7x,'Psi^s(m)')
610   format (i3,4f15.9)
620   format ('  m  k     phi^c(m,k)    phi^c(m,-k)',
     |'     phi^s(m,k)    phi^s(m,-k)')
630   format (2i3,4f15.9)
640   format ('  l  m       q^c(l,m)       q^s(l,m)')
650   format (2i3,2f15.3)
      end
