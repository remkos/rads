* Include file 'GRS80.inc'

      real*8 ae,flat,omfsq,ap,gm,am,ffact,amsq
      parameter (ae=6378137.0d0,flat=1d0/298.257d0)
      parameter (gm=398600.436d9)

* Derived parameters

      parameter (omfsq=(1d0-flat)*(1d0-flat),ap=ae*(1d0-flat))
      parameter (ffact=omfsq-1d0)
      parameter (am=ae*(1d0-flat/3d0-flat*flat/5d0),amsq=am*am)
