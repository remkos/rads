########################################################################
#
#	Makefile for programs handling XGF/XTF/XAF/XXF/XXO data sets
#
# 15-Feb-1998 - Multiplatform version
#
########################################################################

include $(ALTIM)/include/config.mk

########################################################################

LN	= ln -fs

########################################################################

EXEC	= $(BIN)/lxf $(BIN)/ladr $(BIN)/xstat $(BIN)/xtfplot $(BIN)/xxo2xxo \
	$(BIN)/xgf2xgf $(BIN)/xaf2xgf $(BIN)/xxf2xgf \
	$(BIN)/adr2xgf $(BIN)/adr2adr $(BIN)/adr2asc $(BIN)/adrlist \
	$(BIN)/asc2adr $(BIN)/adrplot $(BIN)/earthfit \
	$(BIN)/xttag $(BIN)/xttag-h $(BIN)/xbias $(BIN)/xssb $(BIN)/grd2xgf \
	$(BIN)/xgfsplit

all:	$(EXEC)

clean:
	$(RM) -f $(EXEC) *.o

# Programs that require explicit typing

$(BIN)/lxf lxf:	lxf.o $(RSSUBS)
	$(FC) -o $@ lxf.o $(RSSUBS)
	strip $@

$(BIN)/xgfsplit:	xgfsplit.o
	$(FC) -o $@ xgfsplit.o
	strip $@

# Other, more standard, rules

$(BIN)/ladr ladr: ladr.o $(RSSUBS)
	$(FC) -o $@ ladr.o $(RSSUBS)
	strip $@

$(BIN)/xxo2xxo xxo2xxo:	xxo2xxo.o $(RSSUBS)
	$(FC) -o $@ xxo2xxo.o $(PMPLOT) $(RSSUBS)

$(BIN)/xtfplot:	xtfplot.o $(RSSUBS)
	$(FC) -o $@ xtfplot.o $(PMPLOT) $(RSSUBS)
	strip $@

$(BIN)/xstat xstat:	xstat.o $(RSSUBS) $(GRID)
	$(FC) -o $@ xstat.o $(PMPLOT) $(RSSUBS) $(GRID)
	strip $@

$(BIN)/grd2xgf:	grd2xgf.o xtflimits.o $(RSSUBS) $(GRID)
	$(FC) -o $@ grd2xgf.o xtflimits.o $(RSSUBS) $(GRID)
	strip $@

$(BIN)/xgf2xgf:	xgf2xgf.o xtflimits.o $(RSSUBS) $(GRID)
	$(FC) -o $@ xgf2xgf.o xtflimits.o $(RSSUBS) $(GRID)
	strip $@

$(BIN)/xaf2xgf $(BIN)/xxf2xgf:	$(BIN)/xgf2xgf
	$(LN) xgf2xgf $@

$(BIN)/adr2xgf:	adr2xgf.o normal.o $(NUMLIB) $(RSSUBS) $(GRID)
	$(FC) -o $@ adr2xgf.o normal.o $(NUMLIB) $(RSSUBS) $(GRID)
	strip $@

$(BIN)/adr2adr $(BIN)/adr2asc $(BIN)/adrlist $(BIN)/asc2adr:	$(BIN)/adr2xgf
	$(LN) adr2xgf $@

adrplot $(BIN)/adrplot:	adrplot.o $(RSSUBS)
	$(FC) -o $@ adrplot.o $(PMPLOT) $(RSSUBS)
	strip $@

xttag $(BIN)/xttag:	xttag.o $(RSSUBS) $(GRID)
	$(FC) -o $@ xttag.o $(MATHLIB) $(PMPLOT) $(RSSUBS) $(GRID)
	strip $@

$(BIN)/xbias $(BIN)/xssb $(BIN)/xttag-h:	$(BIN)/xttag
	$(LN) xttag $@

earthfit $(BIN)/earthfit:	earthfit.o $(RSSUBS)
	$(FC) -o $@ earthfit.o $(MATHLIB) $(RSSUBS)
	strip $@

########################################################################
# How to make the format manual

man:	$(DOCDIR)/xformats.pdf
$(DOCDIR)/xformats.pdf:	xformats.tex
	pdflatex xformats
	mv xformats.pdf $@
	rm -f xformats.{dvi,log,toc}
