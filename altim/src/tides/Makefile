########################################################
# Makefile for the tides prediction software
########################################################
ALTIM=/Users/davidtrossman/Documents/Code/rads/altim
HOSTTYPE=x86_64
#HOSTTYPE=i386_gfortran
#NETCDFHOME=/opt/local
NETCDFHOME=/opt/local/lib

include $(ALTIM)/include/config.mk

SHELL	= bash
UTHASH	= $(LIBFES)/third_party/uthash
LIBFES	= fes-2.9.3-Source

INC1	= $(LIBFES)/include/fes.h

OBJS1	= $(LIBFES)/build/src/libfes.a
OBJS2	= libfes.o fes_extra.o
OBJS3	= airtide.o gottide.o festide.o poletide.o hrettide.o etide_ce.o lpetide.o webtide.o tides_aux.o
OBJS4	= ../rssubs/freeunit.o ../rssubs/checkenv.o

PROGS	= etide_test fes2nc fes_test got2nc nc2got testtide

TEXSUBS	= $(OBJS3:.o=.tex)
OBJS	= $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)
MAN		= $(DOCDIR)/tides.pdf
FFLAGS	+=-I$(NETCDFHOME)/include

########################################################

.SUFFIXES: .o .f .F .tex .mod

%.o:	%.mod

###########################################

all:	$(TIDES) $(BINDIR)/uv2cs

$(TIDES): $(OBJS) tides.o
	cp -f $(OBJS1) $(TIDES)
	$(AR) dv $(TIDES) compat.c.o
	$(AR) $(ARFLAGS) $(TIDES) $(OBJS2) $(OBJS3) $(OBJS4)
	$(RANLIB) $@

# Create libfes
$(INC1):
	tar -xJf $(LIBFES).tar.xz

$(OBJS1):	$(INC1)
	(cd $(LIBFES); mkdir -p build; cd build; cmake ..; make)

$(OBJS2):	$(OBJS1)

fes_extra.o:	fes_extra.c $(INC1)
	$(CC) $(CFLAGS) -I$(LIBFES)/include -I$(LIBFES)/src -I$(UTHASH)/src -c -o $@ fes_extra.c

fes_test:	fes_test.f90 $(OBJS1) $(OBJS2)
	$(FC) $(FFLAGS) -o $@ fes_test.f90 $(OBJS1) $(OBJS2) -L$(NETCDFHOME)/lib -lnetcdf

# Tests

test:	do_fes_test do_etide_test do_testtide do_hrettidetest

do_fes_test:	fes_test
	(export FES_DATA=$(LIBFES)/test/data ; fes_test $(LIBFES)/test/fes.ini)

do_etide_test:	etide_test etide_test.in
	time etide_test < etide_test.in

do_testtide:	testtide
	testtide GOT4.10c
	testtide FES2004
	testtide FES2014/standard

do_hrettidetest:	hrettidetest
	hrettidetest $(ALTIM)/data/HRET_v8.1/HRET_v8.1_compressed.nc

# Cleanup

clean:
	@-$(RM) -rf $(OBJS) *.o core *.% $(LIBFES) $(TEXSUBS) $(PROGS)

# Dependencies

flather.o:	flather.inc
testtide.o:	tides.o
lpetide.o gottide.o etice_ce.o:	tides_aux.o

# Module

tides.o:	tides.f90 $(OBJS1) $(OBJS2)
	$(FC) $(FFLAGS) tides.f90 -c -o $@
	mv *.mod $(INCDIR)

# Programs

testtide:	testtide.o $(TIDES)
	$(FC) -o $@ testtide.o $(TIDES) $(NETCDF)
	strip $@

hrettidetest:	hrettidetest.o $(TIDES)
	$(FC) -o $@ hrettidetest.o $(TIDES) $(NETCDF)

ascii2grd_fes:	ascii2grd_fes.o $(TIDES)
	$(FC) -o $@ ascii2grd_fes.o $(TIDES)
	strip $@

ascii2grd_got:	ascii2grd_got.o $(TIDES)
	$(FC) -o $@ ascii2grd_got.o $(TIDES)
	strip $@

fes2nc:	fes2nc.o $(RSSUBS)
	$(FC) -o $@ fes2nc.o $(RSSUBS) $(NETCDF)
	strip $@

got2nc:	got2nc.o $(RSSUBS)
	$(FC) -o $@ got2nc.o $(RSSUBS) $(NETCDF)
	strip $@

nc2fes:	nc2fes.o $(RSSUBS)
	$(FC) -o $@ nc2fes.o $(RSSUBS) $(NETCDF)
	strip $@

nc2got:	nc2got.o $(RSSUBS)
	$(FC) -o $@ nc2got.o $(RSSUBS) $(NETCDF)
	strip $@

grd2ascii_fes:	grd2ascii_fes.o $(GRID)
	$(FC) -o $@ grd2ascii_fes.o $(GRID)
	strip $@

$(BINDIR)/uv2cs:	uv2cs.o
	$(FC) -o $@ uv2cs.o
	strip $@

etide_test:	etide_ce.o etide_test.o tides_aux.o
	$(FC) -o $@ etide_test.o etide_ce.o tides_aux.o
	strip $@

# Manual

man:	$(MAN)

$(MAN):	$(TEXSUBS)
	pdflatex tides
	mv tides.pdf $@
	rm -f tides.{dvi,log}

# Data conversions

air:	$(shell echo $(RADSROOT)/altim/data/airtide/airtide.nc)
got:	$(shell echo $(RADSROOT)/altim/data/GOT{99.2b,00.2,4.7,4.8,4.9,4.10,4.10a,4.10c,4.10c_extrapolated}/{ocean,load}tide.nc)
fes:	$(shell echo $(RADSROOT)/altim/data/FES{95.2.1,99,2004}/{ocean,load}tide.nc)

$(RADSROOT)/altim/data/GOT99.2b/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 ; do \
	zcat $(RADSROOT)/altim/data/GOT99.2b/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/GOT00.2/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 ; do \
	zcat $(RADSROOT)/altim/data/GOT00.2/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/GOT4.7/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 m4 s1 ; do \
	zcat $(RADSROOT)/altim/data/GOT4.7/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/GOT4.8/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 m4 s1 ; do \
	zcat $(RADSROOT)/altim/data/GOT4.8/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/GOT4.9/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 m4 s1 ; do \
	zcat $(RADSROOT)/altim/data/GOT4.9/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/GOT4.10/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 m4 s1 ; do \
	zcat $(RADSROOT)/altim/data/GOT4.10/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/GOT4.10a/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 m4 s1 ; do \
	zcat $(RADSROOT)/altim/data/GOT4.10a/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/GOT4.10c/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 m4 s1 ; do \
	zcat $(RADSROOT)/altim/data/GOT4.10c/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/GOT4.10c_extrapolated/%.nc:	got2nc
	rm -f $@
	(for compo in q1 o1 p1 k1 n2 m2 s2 k2 m4 s1 ; do \
	zcat $(RADSROOT)/altim/data/GOT4.10c_extrapolated/grids_$*/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/airtide/%.nc:	got2nc
	rm -f $@
	(for compo in S1 P1 K1 S2 T2 R2; do \
	zcat $(RADSROOT)/altim/data/airtide/data/$${compo}*.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/DTU10/oceantide.nc:	fes2nc
	rm -f $@
	(for compo in Q1 O1 P1 K1 N2 M2 S2 K2 M4 S1 ; do \
	zcat $(RADSROOT)/altim/data/DTU10/DTU10_TIDEMODEL/OCEAN_TIDE_GRIDS/FES2004_FORMAT/$${compo}DTU10.asc.gz | fes2nc $@ $${compo} ; done )
	(for compo in 2N2 Mf Mm Mtm MSqm; do \
	cat $(RADSROOT)/altim/data/FES2004/data/$${compo}_fes*.asc | fes2nc $@ $${compo} ; done )

$(RADSROOT)/altim/data/DTU10/loadtide.nc:	fes2nc got2nc
	rm -f $@
	(for compo in Q1 O1 P1 K1 N2 M2 S2 K2 M4 S1 ; do \
	zcat $(RADSROOT)/altim/data/DTU10/DTU10_TIDEMODEL/OCEAN{+LOAD,}_TIDE_GRIDS/FES2004_FORMAT/$${compo}DTU10.asc.gz | fes2nc $@ $${compo} ; done )
	(for compo in 2N2; do \
	cat $(RADSROOT)/altim/data/FES2004/data/$${compo}_drfes*.asc | fes2nc $@ $${compo} ; done )
	(for compo in mf mm mtm msqm ; do \
	zcat $(RADSROOT)/altim/data/FES2004/extra/$${compo}load.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/FES2004/oceantide.nc:	fes2nc
	rm -f $@
	(for compo in Q1 O1 K1 2N2 N2 M2 K2 S2 P1 M4 S1 Mf Mm Mtm MSqm; do \
	cat $(RADSROOT)/altim/data/FES2004/data/$${compo}_fes*.asc | fes2nc $@ $${compo} ; done )

$(RADSROOT)/altim/data/FES2004/loadtide.nc:	fes2nc got2nc
	rm -f $@
	(for compo in Q1 O1 K1 2N2 N2 M2 K2 S2 P1; do \
	cat $(RADSROOT)/altim/data/FES2004/data/$${compo}_drfes*.asc | fes2nc $@ $${compo} ; done )
	(for compo in m4 s1 mf mm mtm msqm ; do \
	zcat $(RADSROOT)/altim/data/FES2004/extra/$${compo}load.d.gz | got2nc $@ ; done )

$(RADSROOT)/altim/data/FES2002/oceantide.nc:	fes2nc
	rm -f $@
	(for compo in Q1 O1 K1 2N2 N2 M2 K2 S2 P1; do \
	cat $(RADSROOT)/altim/data/FES2002/data/$${compo}_fes*.asc | fes2nc $@ $${compo} ; done )

$(RADSROOT)/altim/data/FES2002/loadtide.nc:	fes2nc
	rm -f $@
	(for compo in Q1 O1 K1 2N2 N2 M2 K2 S2 P1; do \
	cat $(RADSROOT)/altim/data/FES2002/data/$${compo}_drfes*.asc | fes2nc $@ $${compo} ; done )

$(RADSROOT)/altim/data/FES9%/oceantide.nc:	fes2nc
	rm -f $@
	(for compo in Q1 O1 K1 2N2 N2 M2 K2 S2; do \
	cat $(RADSROOT)/altim/data/FES9$*/data/$${compo}_fes*.asc | fes2nc $@ $${compo} ; done )

$(RADSROOT)/altim/data/FES9%/loadtide.nc:	fes2nc
	rm -f $@
	(for compo in Q1 O1 K1 2N2 N2 M2 K2 S2; do \
	cat $(RADSROOT)/altim/data/FES9$*/data/$${compo}_drfes*.asc | fes2nc $@ $${compo} ; done )
