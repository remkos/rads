* Include file 'covar.inc': parameters etc. for covar.f

      include "GRS80.inc"
      integer mx,nmax
      parameter (mx=10002)
      real*8 ln(0:mx),an(mx),bn(mx)
      common /ccovar/ ln,an,bn,nmax
