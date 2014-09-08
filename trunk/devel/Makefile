#-----------------------------------------------------------------------
# $Id$
#
# Copyright (c) 2011-2014  Remko Scharroo (Altimetrics LLC)
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

include ../config.mk
GRIB := $(shell pkg-config --libs grib_api_f90)
PLM = $(ALTIM)/src/geograph/p_lm.o

# Programs to create
PROGS	= rads_add_ecmwf rads_add_era rads_add_flags rads_add_grid rads_add_ib rads_add_iono rads_add_mog2d \
		rads_add_ncep rads_add_orbit rads_add_refframe rads_add_sla rads_add_ssb rads_add_sst rads_add_surface \
		rads_add_tide rads_add_webtide rads_add_ww3_222 rads_add_ww3_314 \
		rads_fix_c2 rads_fix_reaper rads_fix_sa \
		rads_gen_c2_l1r rads_gen_j2 rads_gen_reaper rads_gen_saral \
		ogdrsplit

# RADS library and linked libraries
LIB     = ../src/librads.a $(ALTIM)/lib/rssubs.a
LDLIBS  = $(NETCDF_LIB) $(GRIB)

# Shell scripts
SCRIPTS = rads_add_common rads_gen_c2_l1r.sh rads_gen_c2_l1r_new.sh \
		rads_gen_saral.sh rads_gen_saral_new.sh rads_sandbox.sh

#-----------------------------------------------------------------------
# What to create?
#-----------------------------------------------------------------------

all:	$(PROGS)

install:	all
	$(INSTALL_DIR) $(bindir) $(datadir)/conf
	$(INSTALL_PROGRAM) $(PROGS) $(bindir)
	$(INSTALL_SCRIPT) $(SCRIPTS) $(bindir)
	$(INSTALL_DATA) *.xml $(datadir)/conf

man:	$(MAN)

#-----------------------------------------------------------------------
# Object dependencies
#-----------------------------------------------------------------------

$(PROGS):	rads_devel.o $(LIB)
rads_add_iono:	$(ALTIM)/lib/iono.a $(ALTIM)/lib/iri2007tec.a
rads_add_ecmwf rads_add_era rads_add_ncep rads_add_tide rads_add_webtide:	$(ALTIM)/lib/tides.a
rads_add_refframe:	$(PLM)

#-----------------------------------------------------------------------
# Test programs
#-----------------------------------------------------------------------

test:

#-----------------------------------------------------------------------
# How to cleanup
#-----------------------------------------------------------------------

clean spotless:
	$(RM) $(PROGS) *.o *.mod
