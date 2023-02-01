* file: flather.inc

* Specify directory name which contains ??.G, ??.H, and tides.inp

      character*80 directory
      parameter (directory='/home/altim/data/flather')

      real*8 pi,twopi,rad
      common /konst2/ pi,twopi,rad

      integer ncomp
      parameter (ncomp=8)
      real*8 freq(ncomp),v(ncomp),u(ncomp),f(ncomp)
      common /compo2/ freq,v,u,f
