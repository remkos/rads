########################################################################
#
#       Makefile for GEODYN Tools
#
# 18-Feb-1998 - Multiplatform version
#
########################################################################

include $(ALTIM)/include/config.mk

########################################################################
# NO ADJUSTMENTS should be needed beyond this point.
########################################################################

EXEC	= $(BIN)/lres $(BIN)/lresconv \
		$(BIN)/g2todr $(BIN)/g2tpoe $(BIN)/g2tsp1 $(BIN)/g2tpon $(BIN)/g2tdut \
		$(BIN)/g2tprc $(BIN)/g2tprl \
		$(BIN)/xxo2asc $(BIN)/asc2gbx \
		$(BIN)/cstg2gbf $(BIN)/gbfsort $(BIN)/gbflist $(BIN)/gbfedit $(BIN)/lgbf \
		$(BIN)/metricedit \
		$(BIN)/gbfgen \
		$(BIN)/ft7merge $(BIN)/giiesumm $(BIN)/giiesum2 \
		$(BIN)/gendrag $(BIN)/genstapos $(BIN)/genaccel9 \
		$(BIN)/autoorbstat $(BIN)/stapos $(BIN)/df2gbf \
		$(BIN)/ran2gbf $(BIN)/dop2met5 \
		$(BIN)/dorisprepro $(BIN)/dorispostpro

all:	$(EXEC)

clean:
		$(RM) -f $(EXEC) *.o

########################################################################

gbfgen $(BIN)/gbfgen:	gbfgen.o $(RSSUBS)
		$(FC) -o $@ gbfgen.o $(RSSUBS)
		strip $@

gbfedit $(BIN)/gbfedit:	gbfedit.o $(RSSUBS)
		$(FC) -o $@ gbfedit.o $(RSSUBS)
		strip $@

metricedit $(BIN)/metricedit:	metricedit.o $(RSSUBS)
		$(FC) -o $@ metricedit.o $(RSSUBS)
		strip $@

gbfedit.o metricedit.o:	gbfedit.inc

df2gbf $(BIN)/df2gbf: df2gbf.o $(RSSUBS)
		$(FC) -o $@ df2gbf.o $(RSSUBS)
		strip $@
df2gbf.o:	df2gbf.f
		$(FC) -O -c df2gbf.f

dorispostpro $(BIN)/dorispostpro:	dorispostpro.o resread.o $(RSSUBS)
		$(FC) -o $@ dorispostpro.o resread.o $(RSSUBS)
		strip $@

dorisprepro $(BIN)/dorisprepro:	dorisprepro.o $(RSSUBS)
		$(FC) -o $@ dorisprepro.o $(RSSUBS)
		strip $@

ran2gbf $(BIN)/ran2gbf: ran2gbf.o
		$(FC) -o $@ ran2gbf.o
		strip $@

dop2met5 $(BIN)/dop2met5: dop2met5.o
		$(FC) -o $@ dop2met5.o
		strip $@

lres $(BIN)/lres:	lres.o resread.o resfit.o resedit.o \
		dopresedt.o resfitedit.o $(RSSUBS)
		$(FC) -o $@ lres.o resread.o resfit.o dopresedt.o resedit.o \
		resfitedit.o $(RSSUBS) $(PMPLOT)
		strip $@
lres.o:		lres.inc
resedit.o:	lres.inc

lresconv $(BIN)/lresconv:	lresconv.o resread.o $(RSSUBS)
		$(FC) -o $@ lresconv.o resread.o $(RSSUBS)
		strip $@

cstg2gbf $(BIN)/cstg2gbf:	catalog.o cstg2gbf.o tdiff.o prerror.o $(RSSUBS)
		$(FC) -o $@ catalog.o cstg2gbf.o tdiff.o prerror.o $(RSSUBS)
		strip $@
catalog.o cstg2gbf.o: cstg2gbf.inc

gbfsort $(BIN)/gbfsort:	gbfsort.o $(RSSUBS)
		$(FC) -o $@ gbfsort.o $(RSSUBS)
		strip $@
gbfsort.o:	gbfsort.inc gbfsort.f

lgbf $(BIN)/lgbf:	lgbf.o
		$(FC) -o $@ lgbf.o
		strip $@

gbflist $(BIN)/gbflist: $(BIN)/gbfsort
		ln -fs gbfsort $@

giiesumm $(BIN)/giiesumm:	giiesumm.o $(RSSUBS)
		$(FC) -o $@ giiesumm.o $(RSSUBS)
		strip $@

giiesum2 $(BIN)/giiesum2:	giiesum2.o $(NUMLIB) $(RSSUBS)
		$(FC) -o $@ giiesum2.o $(NUMLIB) $(RSSUBS)
		strip $@

g2todr $(BIN)/g2todr:	g2todr.o $(RSSUBS)
		$(FC) -o $@ g2todr.o $(RSSUBS)
		strip $@

g2tsp1 $(BIN)/g2tsp1:	g2tsp1.o $(RSSUBS)
		$(FC) -o $@ g2tsp1.o $(RSSUBS)
		strip $@

g2tpoe $(BIN)/g2tpoe:	g2tpoe.o $(RSSUBS)
		$(FC) -o $@ g2tpoe.o $(RSSUBS)
		strip $@

g2tpon $(BIN)/g2tpon:	$(BIN)/g2tpoe
		ln -fs g2tpoe $@

g2tdut $(BIN)/g2tdut:	$(BIN)/g2tpoe
		ln -fs g2tpoe $@

g2tprl $(BIN)/g2tprl:	$(BIN)/g2tpoe
		ln -fs g2tpoe $@

g2tprc $(BIN)/g2tprc:	$(BIN)/g2tpoe
		ln -fs g2tpoe $@

asc2gbx $(BIN)/asc2gbx:	asc2gbx.o $(RSSUBS)
		$(FC) -o $@ asc2gbx.o $(RSSUBS)
		strip $@

xxo2asc $(BIN)/xxo2asc:	xxo2asc.o $(RSSUBS) $(GRID)
		$(FC) -o $@ xxo2asc.o $(RSSUBS) $(GRID)
		strip $@

ft7merge $(BIN)/ft7merge:      ft7merge.o $(RSSUBS)
		$(FC) -o $@ ft7merge.o $(RSSUBS)
		strip $@

gendrag $(BIN)/gendrag:	gendrag.o $(RSSUBS)
		$(FC) -o $@ gendrag.o $(RSSUBS)

genaccel9 $(BIN)/genaccel9:	genaccel9.o $(RSSUBS)
		$(FC) -o $@ genaccel9.o $(RSSUBS)

genstapos $(BIN)/genstapos:	genstapos.o $(RSSUBS)
		$(FC) -o $@ genstapos.o $(RSSUBS)

stapos $(BIN)/stapos:	stapos.o $(RSSUBS)
		$(FC) -o $@ stapos.o $(RSSUBS)

autoorbstat $(BIN)/autoorbstat:	autoorbstat.o $(RSSUBS)
		$(FC) -o $@ autoorbstat.o $(RSSUBS) $(PMPLOT)

########################################################################
# No longer supported
########################################################################

badres $(BIN)/badres:	badres.o resread.o $(RSSUBS)
		$(FC) -o $@ badres.o resread.o $(RSSUBS)
		strip $@

ctlgdays $(BIN)/ctlgdays:	ctlgdays.f
		$(FC) -o $@ ctlgdays.f $(RONLIB)
		strip $@

mrt2gbf $(BIN)/mrt2gbf:	mrt2gbf.f
		$(FC) -o $@ mrt2gbf.f $(RONLIB)
		strip $@
