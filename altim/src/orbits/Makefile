######################################################################
#
#       Makefile for orbit related programs
#
# 14-Feb-1998 - Multiplatform version
#
######################################################################

include $(ALTIM)/include/config.mk

######################################################################

EXEC	= $(BIN)/lodr $(BIN)/odrdiff $(BIN)/odrdif3 $(BIN)/odrhead \
	$(BIN)/odrvalid $(BIN)/esoc2odr $(BIN)/texas2odr $(BIN)/prc2odr \
	$(BIN)/poe2odr $(BIN)/sp12odr $(BIN)/spx2odr $(BIN)/orbite2odr \
	$(BIN)/trackingmap $(BIN)/odr2adr $(BIN)/eof2odr

all:	$(EXEC)

clean:
	$(RM) -f $(EXEC) *.o

######################################################################

$(BIN)/lodr:	lodr.o $(RSSUBS)
	$(FC) -o $@ lodr.o $(RSSUBS)
	strip $@

$(BIN)/odrvalid:	odrvalid.o $(RSSUBS)
	$(FC) -o $@ odrvalid.o $(PMPLOT) $(RSSUBS)
	strip $@

$(BIN)/odrdiff:	odrdiff.o $(RSSUBS) $(SPLIB)
	$(FC) -o $@ odrdiff.o $(PMPLOT) $(RSSUBS) $(SPLIB)
	strip $@

$(BIN)/odrdif3:	odrdif3.o $(RSSUBS) $(SPLIB)
	$(FC) -o $@ odrdif3.o $(PMPLOT) $(RSSUBS) $(SPLIB)
	strip $@

$(BIN)/odrhead:	odrhead.o $(RSSUBS)
	$(FC) -o $@ odrhead.o $(RSSUBS)
	strip $@

$(BIN)/esoc2odr: esoc2odr.o $(RSSUBS)
	$(FC) -o $@ esoc2odr.o $(RSSUBS)
	strip $@

$(BIN)/texas2odr: texas2odr.o $(RSSUBS)
	$(FC) -o $@ texas2odr.o $(RSSUBS)
	strip $@

$(BIN)/poe2odr: poe2odr.o $(RSSUBS)
	$(FC) -o $@ poe2odr.o $(RSSUBS)
	strip $@

$(BIN)/eof2odr: eof2odr.o $(RSSUBS)
	$(FC) -o $@ eof2odr.o $(RSSUBS)
	strip $@

$(BIN)/sp12odr: $(BIN)/spx2odr
	ln -fs spx2odr $@

$(BIN)/spx2odr: spx2odr.o $(RSSUBS) 
	$(FC) -o $@ spx2odr.o $(RSSUBS)
	strip $@

$(BIN)/orbite2odr: orbite2odr.o $(RSSUBS) 
	$(FC) -o $@ orbite2odr.o $(RSSUBS)
	strip $@

$(BIN)/prc2odr: prc2odr.o $(RSSUBS) 
	$(FC) -o $@ prc2odr.o $(RSSUBS)
	strip $@

$(BIN)/odr2adr: odr2adr.o $(RSSUBS) $(GRID)
	$(FC) -o $@ odr2adr.o $(RSSUBS) $(GRID)
	strip $@

$(BIN)/trackingmap: trackingmap.o $(RSSUBS)
	$(FC) -o $@ trackingmap.o $(PMPLOT) $(RSSUBS)
