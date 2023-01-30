########################################################################
#
#	Makefile for grid generating and grid manipulation programs
#
# 15-Feb-1998 - Multiplatform version
#
########################################################################

include $(ALTIM)/include/config.mk
FFLAGS += -I$(ALTIM)/include
BIN = $(ALTIM)/bin

########################################################################

EXEC	= $(BIN)/geogen $(BIN)/gridarea $(BIN)/gridinfo \
		$(BIN)/qgrid $(BIN)/diswgrid $(BIN)/varsgridb \
		$(BIN)/grid2xtra $(BIN)/svd $(BIN)/gridtrend \
		$(BIN)/sstgrid $(BIN)/makessh $(BIN)/lgrid \
		$(BIN)/gridstat $(BIN)/gridamph $(BIN)/hovmuller \
		$(BIN)/dyntopo $(BIN)/gridshift $(BIN)/gridsel \
		$(BIN)/gridinter $(BIN)/grid2nc $(BIN)/xyz2deos \
		$(BIN)/maskinter $(BIN)/svdcmp $(BIN)/basinstat

all:	$(EXEC)

clean:
		$(RM) -f $(EXEC) *.o

########################################################################

dyntopo $(BIN)/dyntopo:	dyntopo.o $(GRID)
		$(FC) -o $@ dyntopo.o $(GRID)
		strip $@

varsgridb $(BIN)/varsgridb:	varsgridb.o $(GRID)
		$(FC) -o $@ varsgridb.o $(GRID)
		strip $@

diswgrid $(BIN)/diswgrid:	diswgrid.o $(GRID) $(RSSUBS)
		$(FC) -o $@ diswgrid.o $(GRID) $(RSSUBS)
		strip $@
diswgrid.o:	diswgrid.f diswgrid.inc

gridsel $(BIN)/gridsel:	gridsel.o $(RSSUBS)
		$(FC) -o $@ gridsel.o $(RSSUBS)
		strip $@

gridarea $(BIN)/gridarea:	gridarea.o
		$(FC) -o $@ gridarea.o
		strip $@

maskinter $(BIN)/maskinter:	maskinter.o $(GRID)
		$(FC) -o $@ maskinter.o $(GRID)
		strip $@

gridinter $(BIN)/gridinter:	gridinter.o
		$(FC) -o $@ gridinter.o $(GRID) $(NETCDF)

xyz2deos $(BIN)/xyz2deos:	xyz2deos.o $(GRID)
		$(FC) -o $@ xyz2deos.o $(GRID)
		strip $@

grid2grid $(BIN)/grid2grid:	grid2grid.o $(GRID)
		$(FC) -o $@ grid2grid.o $(GRID) $(NETCDF)
		strip $@

grid2nc $(BIN)/grid2nc:	grid2nc.o $(GRID) $(RSSUBS)
		$(FC) -o $@ grid2nc.o $(GRID) $(RSSUBS) $(NETCDF)
		strip $@

gridinfo $(BIN)/gridinfo:	gridinfo.o $(GRID)
		$(FC) -o $@ gridinfo.o $(GRID)
		strip $@

gridstat $(BIN)/gridstat:	gridstat.o $(GRID)
		$(FC) -o $@ gridstat.o $(GRID)
		strip $@

gridamph $(BIN)/gridamph:	gridamph.o $(GRID)
		$(FC) -o $@ gridamph.o $(GRID)
		strip $@

grid2xtra $(BIN)/grid2xtra:	grid2xtra.o $(GRID)
		$(FC) -o $@ grid2xtra.o $(GRID)
		strip $@

gridshift $(BIN)/gridshift:	gridshift.o $(GRID)
		$(FC) -o $@ gridshift.o $(GRID)
		strip $@

lgrid $(BIN)/lgrid:	lgrid.o $(GRID)
		$(FC) -o $@ lgrid.o $(GRID)
		strip $@

geogen $(BIN)/geogen:	geogen.o legpol.o $(GRID) $(RSSUBS)
		$(FC) -o $@ geogen.o legpol.o $(GRID) $(RSSUBS)
		strip $@

sstgrid $(BIN)/sstgrid:	sstgrid.o legpol.o $(GRID) $(RSSUBS)
		$(FC) -o $@ sstgrid.o legpol.o $(GRID) $(RSSUBS)
		strip $@

makessh $(BIN)/makessh:	makessh.o $(GRID) $(RSSUBS)
		$(FC) -o $@ makessh.o $(GRID) $(RSSUBS)
		strip $@

hovmuller $(BIN)/hovmuller:	hovmuller.o $(GRID) $(RSSUBS)
		$(FC) -o $@ hovmuller.o $(GRID) $(RSSUBS)
		strip $@

gridtrend $(BIN)/gridtrend:	gridtrend.o $(GRID) $(RSSUBS)
		$(FC) -o $@ gridtrend.o $(GRID) $(RSSUBS)
		strip $@

svd $(BIN)/svd:		svd.o dsvdcmp.o dpythag.o $(GRID) $(RSSUBS)
		$(FC) -o $@ svd.o dsvdcmp.o dpythag.o $(GRID) $(RSSUBS)
		strip $@

svdcmp $(BIN)/svdcmp:		svdcmp.o dsvdcmp.o dpythag.o $(GRID) $(RSSUBS)
		$(FC) -o $@ svdcmp.o dsvdcmp.o dpythag.o $(GRID) $(RSSUBS)
		strip $@

basinstat $(BIN)/basinstat:	basinstat.o $(GRID) $(RSSUBS)
		$(FC) -o $@ basinstat.o $(GRID) $(RSSUBS) $(NETCDF)

# qgrid may need F90
include $(HOSTTYPE).mk

qgrid1.o:	qgrid.inc

basinstat.o gridinter.o grid2nc.o maskinter.o:	$(ALTIM)/include/gridlib.mod
