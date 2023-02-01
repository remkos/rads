# This file provides some of the variables used in Makefiles
# for the Altimetry software at DUT/DEOS
#
# Version 0301.0 - Remko Scharroo
#######################################################################
# Object locations: Library, Manuals, Binaries
export ALTIM=/Users/davidtrossman/Documents/Code/rads/altim
export HOSTTYPE=x86_64
export NETCDFHOME=/opt/local

BINDIR	= $(ALTIM)/bin
DOCDIR	= $(ALTIM)/doc
INCDIR	= $(ALTIM)/include
LIBDIR	= $(ALTIM)/lib
SRCDIR	= $(ALTIM)/src
BIN	= $(BINDIR)

# Libraries
GDRSUBS	= $(LIBDIR)/gdrsubs.a
GRID	= $(LIBDIR)/grid.a
NUMLIB	= $(LIBDIR)/numlib.a
PGPLOT	= -L$(LIBDIR) -lpgplot -lX11
PMPLOT	= -L$(LIBDIR) -lpmplot -lX11
PVSUBS	= $(LIBDIR)/pvsubs.a
RECIPES	= $(LIBDIR)/recipes.a
RSSUBS	= $(LIBDIR)/rssubs.a
SPLIB	= $(LIBDIR)/splib.a
TIDES	= $(LIBDIR)/tides.a
IONO	= $(LIBDIR)/iono.a
IRI95	= $(LIBDIR)/iri95tec.a
IRI2007	= $(LIBDIR)/iri2007tec.a
RADS	= $(RADSROOT)/lib/rads.a
NETCDF	?= $(if $(NETCDFHOME),-L$(NETCDFHOME)/lib -lnetcdf)

# .F to .o mechanism
.SUFFIXES:	.F
.F.o:
	$(FC) $(FFLAGS) $< -c

# .f90 rules
.SUFFIXES:	.f90
.f90.o:
	$(COMPILE.f) $(OUTPUT_OPTION) $<
.f90:
	$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

# .tex and .dvi mechanisms
.SUFFIXES:	.tex .dvi
.c.tex:
	maketex $< > $@
.f.tex:
	maketex $< > $@
.F.tex:
	maketex $< > $@
.f90.tex:
	maketex $< > $@
.tex.dvi:
	latex $<

# C preprocessor and other executables
CPP	= /lib/cpp -P $(IFLAGS) $(CPPFLAGS)
RANLIB	= ranlib
AR	= ar

# First some system dependencies
include $(INCDIR)/$(HOSTTYPE)/sysdep.mk

# Fortran/C Compiler flags
IFLAGS  = -I. -I$(INCDIR) -I$(INCDIR)/$(HOSTTYPE) $(if $(NETCDFHOME),-I$(NETCDFHOME)/include)
FFLAGS += $(IFLAGS) -O2
CFLAGS += $(IFLAGS) -O2
