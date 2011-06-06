########################################################################
#
#       Makefile for RADS4 library and utilities
#
# $Id: $
########################################################################

include config.mk

%.o:	%.mod

OBJS1	= rads_netcdf.o rads_time.o rads_math.o rads_misc.o rads_grid.o xmlparse.o
OBJS2	= rads.o rads3.o
OBJS	= $(OBJS1) $(OBJS2)

LIB	= librads.a

PROGS_O = rads2asc.o rads2adr.o rads2nc.o rads2grd.o radsstat.o radscolin.o radsxogen.o radsxoconv.o \
	rads_test.o rads_nml2xml.o rads_sample.o rads3test.o
PROGS_L = rads2xadr rads2xgf
PROGS_L4= rads2asc4 radscolin4 radsstat4 rads2grd4
PROGS = $(PROGS_O:.o=) $(PROGS_L) $(PROGS_L4)

INCS	= config.inc

MODULES = typesizes.mod $(OBJS:.o=.mod)

########################################################################
# What to create?
########################################################################

all:	$(LIB) $(PROGS)

install:	all
	$(INSTALL) -d $(LIBDIR)
	$(INSTALL) -m 644 $(LIB) $(LIBDIR)
	$(INSTALL) -d $(INCDIR)
	$(INSTALL) -m 644 ../include/*.h ../include/*.inc ../include/*.mk $(INCDIR)
	$(INSTALL) -m 644 *.mod $(INCDIR)
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 755 $(PROGS) $(BINDIR)
	( cd $(BINDIR) ; ln -fs rads2adr rads2xadr )
	( cd $(BINDIR) ; ln -fs rads2adr rads2xgf )

man:	$(MAN)

########################################################################
# Object dependencies
########################################################################

$(PROGS_O): rads.mod rads_misc.mod rads_time.mod rads_netcdf.mod
rads_nml2xml.o:	xmlparse.mod

nag.o:  nag.f
	$(FC) -O -c nag.f

$(OBJS):	typesizes.mod
rads.o:		$(INCS) $(OBJS1)
rads3.o rads_test.o:	rads.o

########################################################################
# Create config.inc when config.mk is updated
########################################################################
config.inc:	config.mk
	printf "character(80) :: radsdataroot = '%s', radsuserroot = '%s'\n" $(datarootdir) $(HOME) > $@

########################################################################
# Create Fortran 90 interface if needed
########################################################################

typesizes.mod:	typeSizes.f90
	$(FC) $(FFLAGS) -c typeSizes.f90
	rm -f typeSizes.o

########################################################################
# How to create the library
########################################################################

$(LIB):	$(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@

########################################################################
# Test programs
########################################################################
test:	rads_test
	rads_test debug=3 sat=e1/a ymd=19910924,19910925 lat=-50,30 lon=-90,120 \
	sel=time,lat,lon,sla,12%,tide_ocean%,pass,drange > $(TMPDIR)/rads_test.out
	diff rads_test.out $(TMPDIR)/rads_test.out
rads_test:	$(LIB) rads_test.o
	$(FC) -o $@ rads_test.o $(LIB) $(NETCDF)
rads_sample:	$(LIB) rads_sample.o
	$(FC) -o $@ rads_sample.o $(LIB) $(NETCDF)
rads_nml2xml:	$(LIB) rads_nml2xml.o
	$(FC) -o $@ rads_nml2xml.o $(LIB)
rads3test:	$(LIB) rads3test.o
	$(FC) -o $@ rads3test.o $(LIB) $(NETCDF)

########################################################################
#       Definition of executable generation
########################################################################

radsconfig.sh:	radsconfig.sh.in
	@echo "Error: Need to run configure in src directory first"
	@exit 1

cyclemean:	cyclemean.o $(LIB)
	$(FC) -o $@ cyclemean.o $(LIB) $(NETCDF)

rads2asc: rads2asc.o $(LIB)
	$(FC) -o $@ rads2asc.o $(LIB) $(NETCDF)

rads2adr: rads2adr.o $(LIB)
	$(FC) -o $@ rads2adr.o $(LIB) $(NETCDF)

rads2nc: rads2nc.o $(LIB)
	$(FC) -o $@ rads2nc.o $(LIB) $(NETCDF)

$(PROGS_L): rads2adr
	ln -fs rads2adr $@

$(PROGS_L4):	%4:	%
	ln -fs $* $@

rads2grd:	rads2grd.o $(LIB)
	$(FC) -o $@ rads2grd.o $(LIB) $(NETCDF)

radsstat:	radsstat.o $(LIB)
	$(FC) -o $@ radsstat.o $(LIB) $(NETCDF)

radscolin:	radscolin.o $(LIB)
	$(FC) -o $@ radscolin.o $(LIB) $(NETCDF)

radscont:	radscont.o $(LIB)
	$(FC) -o $@ radscont.o $(LIB) $(NETCDF)

radsxogen:	radsxogen.o iqsort.o nag.o $(LIB)
	$(FC) -o $@ radsxogen.o iqsort.o nag.o $(LIB) $(NETCDF)

radsxoconv:	radsxoconv.o $(LIB)
	$(FC) -o $@ radsxoconv.o $(LIB) $(NETCDF)
