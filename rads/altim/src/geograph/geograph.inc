*
* file: geograph.inc
*
      real*8 pi,rad
      parameter (pi=3.14159265358979d0,rad=pi/180d0)

      real*8 gm,ae,flat,a0,incl,n0,wmdot,ogdot
      common /cgeogr1/ gm,ae,flat,a0,incl,n0,wmdot,ogdot

      real*8 deadband
      logical test
      common /cgeogr2/ deadband,test

      integer ndeg,narray,npar,ncoef
      parameter (ndeg=120)
      parameter (narray=(2*ndeg+3)*(ndeg+2)*(ndeg+1)/6-1)
      parameter (npar=ndeg*(ndeg+2),ncoef=(ndeg+1)*(ndeg+2)/2)
      integer pnt(0:ndeg)
      real*8 cosi(0:ndeg),sini(0:2*ndeg),bin(0:ndeg,0:ndeg),
     |       qc(ncoef),qs(ncoef),dlmp(0:narray),
     |       cosml(0:ndeg),sinml(0:ndeg)
      common /cgeogr3/ cosi,sini,qc,qs,dlmp,cosml,sinml,bin,pnt

      real*8 cs(npar)
      integer ics(npar),iord(npar),ideg(npar),h(0:ndeg),lmax,ipmax
      common /cgeogr4/ cs,ics,iord,ideg,h,lmax,ipmax

      real*4 covar(npar*(npar+1)/2)
      real*4 solve(npar),sig(npar),rvect(npar)
      integer g(npar)
      logical docovar
      common /cgeogr5/covar,solve,sig,rvect,g,docovar
