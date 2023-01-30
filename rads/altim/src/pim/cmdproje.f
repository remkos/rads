      subroutine cmdproje
      include "pim.inc"
      real trspar(2)/-100.,-100./

* Set projection.
* If you specify a projection type 0, NO GEOGRAPHICAL PROJECTION is used,
* but rather a standard PGSWIN/PGBOX environment.

      project=1
      if (popcmd('PROJE',argum)) then
	 call pimopt('tsp=',argum,trspar(1),trspar(2),dum,dum)
	 read (argum,*,iostat=ios) project,trspar
      endif
      call pmdef(project,0.0,trspar(1),trspar(2))

      end
