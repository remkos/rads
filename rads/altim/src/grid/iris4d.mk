########################################################################
#
#	Makefile for qgrid
#
########################################################################

qgrid.o:	qgrid.F qgrid.inc
		$(F90) $(FFLAGS) -Df90 qgrid.F -c

qgrid $(BIN)/qgrid:	qgrid.o qgrid1.o ../xover/xtflimits.o $(GRID) $(RSSUBS)
		$(F90) -o $@ qgrid.o qgrid1.o ../xover/xtflimits.o $(GRID) $(RSSUBS)
