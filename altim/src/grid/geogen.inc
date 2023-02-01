
c-----------------------------------------------------------------------
c parameter definitions file
c-----------------------------------------------------------------------
      integer absmax,cmax,rootsize,coefsize,foursize
      parameter (pi       = 3.141592653589793d0)
      parameter (rootsize = 2000)      ! max. for root computation
      parameter (absmax   = 360)        ! max. array dimension
      parameter (cmax     = 360)       ! max.degree of pot.model
      parameter (coefsize = (cmax+1)*(cmax+2)/2) ! max. nr. of pot. coeff.

      parameter (lucoef   = 1)         ! pot. coeff. model input unitnr.
      parameter (lugra    = 2)         ! gradient info output unitnr.
      parameter (in       = 5)         ! user input. 
      parameter (out      = 6)         ! program output.

      parameter (mmin     =  0)

c-----------------------------------------------------------------------
c variables and arrays definition file
c--------------------------------------------------------------------
      real*8  cnm(coefsize),snm(coefsize)
      real*8  gm,ae,root(0:rootsize),height
      integer offset(0:coefsize)

c-----------------------------------------------------------------------
c common block definitions file
c-----------------------------------------------------------------------
      common /ptcoef/ cnm,snm,offset
      common /rtdata/ root
      common /fsdata/ gm,ae,height
      common /input/ nmax,lmin
