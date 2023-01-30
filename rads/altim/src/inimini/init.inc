* init.inc

* edtmlt	residual edit multiplier
* edttrk	track edit multiplier
* perc		percentage convergence criterion
* maxvar	maximum sea surface variability
* maxiter	maximum number of iterations
* smcells	number of cells in the smoothing process of the rms grid

      real    edtmlt,maxvar,perc,version,edttrk
      integer maxiter,smcells

      common /misc/  edtmlt,perc,edttrk,maxvar,version,maxiter,smcells

      real    rad,pi,sqrt_two,murad

      common /math/  pi,rad,sqrt_two,murad

      integer wrap(-32768:32767)

      common /cwrap/ wrap
