########################################################################
#
#       Makefile for GEOGRAPH
#
# 07-Apr-1998 - Multiplatform version
#
########################################################################
include $(ALTIM)/include/config.mk

########################################################################
# NO ADJUSTMENTS should be needed beyond this point.
########################################################################
.SUFFIXES:	dd.o .F

GEOSUBS = f_lmp.o p_lm.o d_lmp.o q_lm.o gravrd.o geograph.o geogrcmp.o mangle.o

EXEC	= $(BIN)/geospec $(BIN)/geomean $(BIN)/geogrid $(BIN)/geoint \
	$(BIN)/geox $(BIN)/geox1 $(BIN)/geox2 $(BIN)/lcovar $(BIN)/geogen3 \
	$(BIN)/xxosort

all:	$(EXEC)
extra:	gravdiff gravcorr readcov rewrcov rewrnor ftest geogrtst
clean:
	$(RM) -f $(EXEC) *.o

$(GEOSUBS) gravdiff.o geogen3.o geogrid.o geogrtst.o geoint.o geomean.o \
geospec.o geox1.o geox2.o geox.o gravdiff.o gravrd.o: geograph.inc

.Fdd.o:
	$(FC) $(FFLAGS) -DDIFF -c -o $*dd.o $*.F

.F.f:
	$(CPP) < $*.F > $*.f

$(GEOSUBS) geox.o geospecc.o geospec.o geoint.o geomean.o geogrid.o: geograph.inc

geosubs.a:	$(GEOSUBS)
		-rm -f $@
		$(AR) $(ARFLAGS) $@ $(GEOSUBS)
		$(RANLIB) $@

ftest:		ftest.o f_lmp.o p_lm.o mangle.o
		$(FC) -o $@ ftest.o f_lmp.o p_lm.o mangle.o

ftestdd:	ftestdd.o f_lmpdd.o p_lmdd.o mangle.o
		$(FC) -o $@ ftestdd.o f_lmpdd.o p_lmdd.o mangle.o

$(BIN)/geomean:	geomean.o geosubs.a $(RSSUBS)
		$(FC) -o $@ geomean.o geosubs.a $(RSSUBS)

$(BIN)/geoint:	geoint.o geosubs.a $(RSSUBS)
		$(FC) -o $@ geoint.o geosubs.a $(RSSUBS)

$(BIN)/geospec:	geospec.o geosubs.a $(RSSUBS)
		$(FC) -o $@ geospec.o geosubs.a $(PMPLOT) $(RSSUBS)

$(BIN)/geox:	geox.o gridrms.o geosubs.a $(GRID) $(RSSUBS)
		$(FC) -o $@ geox.o gridrms.o geosubs.a $(GRID) $(RSSUBS) $(MATHLIB)

$(BIN)/geox1:	geox1.o geosubs.a $(RSSUBS)
		$(FC) -o $@ geox1.o geosubs.a $(RSSUBS)

$(BIN)/geox2:	geox2.o geosubs.a $(GRID) $(RSSUBS)
		$(FC) -o $@ geox2.o geosubs.a $(GRID) $(RSSUBS) $(MATHLIB)

$(BIN)/geogrid:	geogrid.o gridrms.o geosubs.a $(GRID) $(RSSUBS)
		$(FC) -o $@ geogrid.o gridrms.o geosubs.a $(GRID) $(RSSUBS)

$(BIN)/lcovar:	lcovar.o
		$(FC) -o $@ lcovar.o

$(BIN)/geogen3:	geogen3.o geosubs.a $(GRID) $(RSSUBS)
		$(FC) -o $@ geogen3.o geosubs.a $(GRID) $(RSSUBS)

$(BIN)/xxosort:	xxosort.o $(GRID) $(RSSUBS)
		$(FC) -o $@ xxosort.o $(GRID) $(RSSUBS) $(MATHLIB)

geogrtst:	geogrtst.o d_lmp.o p_lm.o f_lmp.o q_lm.o geograph.o mangle.o
		$(FC) -o $@ geogrtst.o d_lmp.o p_lm.o f_lmp.o q_lm.o geograph.o mangle.o

gravcorr:	gravcorr.o $(GRID) $(RSSUBS)
		$(FC) -o $@ gravcorr.o $(GRID) $(RSSUBS)

gravdiff:	gravdiff.o geosubs.a
		$(FC) -o $@ gravdiff.o geosubs.a

rewrnor:	rewrnor.o
		$(FC) -o $@ rewrnor.o

rewrcov:	rewrcov.o
		$(FC) -o $@ rewrcov.o

readcov:	readcov.o
		$(FC) -o $@ readcov.o

test.o:		geograph.inc
test:		test.o gravrd.o
		$(FC) -o $@ test.o gravrd.o
