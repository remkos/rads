########################################################################
#
#	Makefile for qgrid
#
########################################################################

qgrid.o:	qgrid.F qgrid.inc
		$(FC) $(FFLAGS) $*.F -c

qgrid $(BIN)/qgrid:	qgrid.o qgrid1.o ../xover/xtflimits.o $(GRID) $(RSSUBS)
		$(FC) -o $@ qgrid.o qgrid1.o ../xover/xtflimits.o $(GRID) $(RSSUBS)
