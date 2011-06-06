#######################################################################
# This file provides some of the variables used in Makefiles
# for the RADS software
#
# Version 0304.0 - Remko Scharroo
#######################################################################
# Object locations: Library, Manuals, Binaries
prefix = /rads
exec_prefix = ${prefix}
datarootdir = /rads/data
BINDIR = ${exec_prefix}/bin
DOCDIR = ${datarootdir}/doc/${PACKAGE_TARNAME}
INCDIR = ${prefix}/include
LIBDIR = ${exec_prefix}/lib
RADSROOT = $(prefix)
RADSDATAROOT = $(datarootdir)
RADS   = $(LIBDIR)/librads.a
INSTALL = /usr/bin/install -c

# Include and load option for NetCDF
NETCDF_INC = /sw/lib/netcdf-gfortran/include
NETCDF = -L/sw/lib -L/sw/lib/netcdf-gfortran/lib -lnetcdff -lnetcdf

# Other libraries for tools only
GDRSUBS= $(ALTIM)/lib/gdrsubs.a
GRID   = $(ALTIM)/lib/grid.a
NUMLIB = $(ALTIM)/lib/numlib.a
PGPLOT = -L$(ALTIM)/lib -lpgplot -L/usr/X11/lib -lX11 -lpng
PMPLOT = -L$(ALTIM)/lib -lpmplot -L/usr/X11/lib -lX11 -lpng
PVSUBS = $(ALTIM)/lib/pvsubs.a
RECIPES= $(ALTIM)/lib/recipes.a
RSSUBS = $(ALTIM)/lib/rssubs.a
SPLIB  = $(ALTIM)/lib/splib.a
TIDES  = $(ALTIM)/lib/tides.a
IONO   = $(ALTIM)/lib/iono.a
IRI95  = $(ALTIM)/lib/iri95tec.a
IRI2007= $(ALTIM)/lib/iri2007tec.a

# Modules
MODULES= $(INCDIR)/rads.mod $(INCDIR)/typesizes.mod

# .F to .o rules
.SUFFIXES:	.F
.F.o:
	$(FC) $(FFLAGS) $< -c

# .tex and .dvi rules
.SUFFIXES:	.tex .dvi
.c.tex:
	maketex $< > $@
.f.tex:
	maketex $< > $@
.F.tex:
	maketex $< > $@
.tex.dvi:
	latex $<

# .f90 rules
.SUFFIXES:	.f90
.f90.o:
	$(COMPILE.f) $(OUTPUT_OPTION) $<
.f90:
	$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

# C preprocessor and other executables
RANLIB	= ranlib
AR	= ar
STRIP	= strip

# First some system dependencies
FC	= gfortran
FFLAGS  = -g -O2 -ffixed-line-length-none -Wall -Wtabs -fimplicit-none -I../src -I$(INCDIR) -I$(NETCDF_INC)
#CFLAGS  = -g -O2 -Wall -Wimplicit -I../src -I$(INCDIR) -I$(NETCDF_INC)
MATHLIB = -llapack -lblas
