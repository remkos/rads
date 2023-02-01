########################################################################
#
#       Makefile for INIMINI, INIEDIT, INISURF
#
# 15-Feb-1998 - Multiplatform version
#
########################################################################

include $(ALTIM)/include/config.mk
CPPFLAGS	=

########################################################################

SUBS	= fill_mat.o get_size.o grid.o iniinit.o inimini.o inistat.o \
	mpnt.o rdxaf.o rdxtf.o rdxxf.o solve_mat.o sort.o \
	surface.o update.o vpnt.o wrstat.o wrxaf.o wrxtf.o wrxxf.o xovers.o
INCS	= buffer.inc data.inc grid.inc init.inc satcat.inc stat.inc
EXEC	= $(BIN)/inimini $(BIN)/inisurf $(BIN)/iniedit $(BIN)/xtfrewr

########################################################################

all:	$(EXEC)

$(SUBS) xtfrewr.o:	$(INCS)

$(BIN)/inimini inimini: $(SUBS) $(GRID) $(RSSUBS)
	$(FC) -o $@ $(SUBS) $(GRID) $(MATHLIB) $(RSSUBS)
	strip $@

$(BIN)/inisurf $(BIN)/iniedit: $(BIN)/inimini
	rm -f $@
	ln -fs inimini $@

$(BIN)/xtfrewr: xtfrewr.o $(RSSUBS)
	$(FC) -o $@ xtfrewr.o $(RSSUBS)
	strip $@

clean:
	$(RM) -f $(EXEC) *.o core *.bak *~
