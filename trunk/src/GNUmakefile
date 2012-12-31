#-----------------------------------------------------------------------
# $Id$
#
# Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
# See LICENSE.TXT file for copying and redistribution conditions.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#-----------------------------------------------------------------------

include config.mk

OBJS1	= rads_netcdf.o rads_time.o rads_math.o rads_misc.o rads_grid.o rads_geo.o xmlparse.o
OBJS2	= rads.o
OBJS3	= rads3.o
OBJS	= $(OBJS1) $(OBJS2) $(OBJS3)

MODS1	= $(OBJS1:.o=.mod)
MODS2	= $(OBJS1:.o=.mod)
MODS3	= $(OBJS1:.o=.mod)
MODS	= $(OBJS:.o=.mod)

LIB	= librads.a

PROGS_O = rads2asc.o rads2adr.o rads2nc.o rads2grd.o radsstat.o radscolin.o radsxogen.o radsxoconv.o \
	radsxolist.o rads_test.o radsreconfig.o rads_sample.o rads3test.o radsvar.o
PROGS_L1 = rads2xadr rads2xgf
PROGS_L2 = radsxostat
PROGS_L4 = rads2asc4 radscolin4 radsstat4 rads2grd4
PROGS_L  = $(PROGS_L1) $(PROGS_L2) $(PROGS_L4)
PROGS = $(PROGS_O:.o=)

INCS	= config.inc

#-----------------------------------------------------------------------
# What to create?
#-----------------------------------------------------------------------

all:	$(LIB) $(PROGS) $(PROGS_L)

install:	all
	$(INSTALL) -d $(LIBDIR)
	$(INSTALL) -m 644 $(LIB) $(LIBDIR)
	$(INSTALL) -d $(INCDIR)
	$(INSTALL) -m 644 *.inc *.mod $(INCDIR)
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 755 $(PROGS) $(BINDIR)
	for p in $(PROGS_L1); do ( cd $(BINDIR) ; ln -fs rads2adr $$p ) ; done
	for p in $(PROGS_L2); do ( cd $(BINDIR) ; ln -fs radsxolist $$p ) ; done
	for p in $(PROGS_L4); do ( cd $(BINDIR) ; ln -fs $${p/4/} $$p ) ; done

man:	$(MAN)

#-----------------------------------------------------------------------
# Object dependencies
#-----------------------------------------------------------------------

$(OBJS2):	$(MODS1) $(INCS)
$(OBJS3):	$(MODS2)
$(PROGS_O):	$(MODS)
$(PROGS):	$(LIB)
radsxogen.o:	spline.mod

#-----------------------------------------------------------------------
# Create config.inc when config.mk is updated
#-----------------------------------------------------------------------
config.inc:	config.mk
	printf "character(len=rads_naml) :: radsdataroot = '%s', radsuserroot = '%s'\n" $(datarootdir) $(HOME) > $@

#-----------------------------------------------------------------------
# How to create the library
#-----------------------------------------------------------------------

$(LIB):	$(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@

#-----------------------------------------------------------------------
# Test programs
#-----------------------------------------------------------------------
test:	rads_test
	rads_test debug=3 sat=e1/a ymd=19910924,19910925 lat=-50,30 lon=-90,120 \
	sel=time,lat,lon,sla,12_,tide_ocean_,pass,drange | diff rads_test.out -

rads_test:	rads_test.o
	$(FC) -o $@ rads_test.o $(LIB) $(NETCDF)

rads_sample:	rads_sample.o
	$(FC) -o $@ rads_sample.o $(LIB) $(NETCDF)

radsreconfig:	radsreconfig.o
	$(FC) -o $@ radsreconfig.o $(LIB) $(NETCDF)

rads3test:	rads3test.o
	$(FC) -o $@ rads3test.o $(LIB) $(NETCDF)

#-----------------------------------------------------------------------
#       Definition of executable generation
#-----------------------------------------------------------------------

radsconfig.sh:	radsconfig.sh.in
	@echo "Error: Need to run configure in src directory first"
	@exit 1

rads2asc:	rads2asc.o
	$(FC) -o $@ rads2asc.o $(LIB) $(NETCDF)

rads2adr:	rads2adr.o
	$(FC) -o $@ rads2adr.o $(LIB) $(NETCDF)

rads2nc:	rads2nc.o
	$(FC) -o $@ rads2nc.o $(LIB) $(NETCDF)

rads2grd:	rads2grd.o
	$(FC) -o $@ rads2grd.o $(LIB) $(NETCDF)

radsstat:	radsstat.o
	$(FC) -o $@ radsstat.o $(LIB) $(NETCDF)

radscolin:	radscolin.o
	$(FC) -o $@ radscolin.o $(LIB) $(NETCDF)

radsxogen:	radsxogen.o spline.o
	$(FC) -o $@ radsxogen.o spline.o $(LIB) $(NETCDF)

radsxoconv:	radsxoconv.o
	$(FC) -o $@ radsxoconv.o $(LIB) $(NETCDF)

radsxolist:	radsxolist.o
	$(FC) -o $@ radsxolist.o $(LIB) $(NETCDF)

radsvar:	radsvar.o
	$(FC) -o $@ radsvar.o $(LIB) $(NETCDF)

$(PROGS_L1):	rads2adr
	ln -fs rads2adr $@

$(PROGS_L2):	radsxolist
	ln -fs radsxolist $@

$(PROGS_L4):	%4:	%
	ln -fs $* $@
