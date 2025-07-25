#-----------------------------------------------------------------------
# Copyright (c) 2011-2025  Remko Scharroo
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

# External objects
GRIB_LIB := $(shell pkg-config --libs eccodes_f90 2>/dev/null || pkg-config --libs grib_api_f90 2>/dev/null)
GRIB_INC := $(shell pkg-config --cflags eccodes_f90 2>/dev/null || pkg-config --cflags grib_api_f90 2>/dev/null)
NETCDF_C_LIB	:= $(shell nc-config --libs)
PLM = $(ALTIM)/src/geograph/p_lm.o

# Programs to create
PROGS	= rads_add_dac rads_add_dual rads_add_ecmwf rads_add_era rads_add_era5 rads_add_flags \
		rads_add_gpd rads_add_grid rads_add_ib rads_add_iono rads_add_ncep rads_add_orbit \
		rads_add_refframe rads_add_seaice rads_add_sla rads_add_ssb rads_add_sst rads_add_surface \
		rads_add_tide rads_add_webtide rads_add_mfwam rads_add_ww3_222 rads_add_ww3_314 \
		rads_fix_c2 rads_fix_eqtime rads_fix_j1 rads_fix_jason rads_fix_reaper rads_fix_s3 rads_fix_s6 \
		rads_gen_c2_op rads_gen_h2 rads_gen_jason rads_gen_jason_gdrf rads_gen_reaper rads_gen_s3 rads_gen_s6 \
		rads_gen_saral rads_gen_swot rads_gen_tp_gdrf make_orf ogdrsplit rads_pre_sort_passes rads_copy_var

# Includes
INCLUDES = $(MOD_FLAG). $(MOD_FLAG)../src $(MOD_FLAG)$(ALTIM)/include $(NETCDF_INC) $(GRIB_INC)

# RADS library and linked libraries
LIB     = ../src/librads.a $(ALTIM)/lib/rssubs.a
LDLIBS  = $(NETCDF_LIB) $(GRIB_LIB)

# Shell or perl scripts
SCRIPTS = rads_add_common rads_add_common_gdrg.sh rads_add_era5.sh \
		rads_gen_3a.sh rads_gen_3a_new.sh rads_gen_3a_new_calval.sh rads_gen_3a_calval.sh \
		rads_gen_3b.sh rads_gen_3b_new.sh rads_gen_3b_new_calval.sh rads_gen_3b_calval.sh \
		rads_gen_6a_calval.sh rads_gen_6b_calval.sh \
		rads_gen_c2_op.sh rads_gen_h2.sh \
		rads_gen_j1.sh rads_gen_j2.sh \
		rads_gen_j3.sh rads_gen_j3_calval.sh rads_gen_j3_new.sh rads_gen_j3_new_calval.sh \
		rads_gen_tp_gdrf.sh rads_gen_saral.sh rads_gen_saral_new.sh \
		rads_gen_sw.sh rads_gen_sw_new.sh rads_gen_sw_new_calval.sh rads_gen_sw_calval.sh \
		rads_purge_older.sh rads_sandbox.sh hardlink.sh \
		radscyclesindex

#-----------------------------------------------------------------------
# What to create?
#-----------------------------------------------------------------------

all:	$(PROGS)

install:	all
	$(INSTALL_DIR) $(bindir) $(datadir)/conf
	$(INSTALL_PROGRAM) $(PROGS) $(bindir)
	$(INSTALL_SCRIPT) $(SCRIPTS) $(bindir)
	$(INSTALL_DATA) *.xml $(datadir)/conf
	(cd $(bindir) ; $(LN_S) rads_add_dac rads_add_mog2d)

man:	$(MAN)

#-----------------------------------------------------------------------
# Object dependencies
#-----------------------------------------------------------------------

rads_gen_h2 rads_gen_jason rads_gen_jason_gdrf rads_gen_saral rads_gen_s3 rads_gen_s6 rads_gen_swot rads_gen_tp_gdrf rads_gen_c2_op:	rads_devel_netcdf.o rads_gen.o
rads_gen_reaper:	rads_gen.o
rads_fix_s6: rads_devel.o rads_devel_netcdf.o
rads_add_iono:	$(ALTIM)/lib/iono.a
rads_add_ecmwf rads_add_era rads_add_ncep rads_add_tide rads_add_webtide:	$(ALTIM)/lib/tides.a
rads_add_tide:	LDLIBS = $(NETCDF_LIB) $(NETCDF_C_LIB) $(GRIB_LIB)
rads_add_refframe:	$(PLM)
rads_pre_sort_passes ogdrsplit:	rads_devel_misc.o
rads_devel.o:	../src/rads.mod ../src/rads_misc.mod
rads_devel_misc.o:	../src/rads_misc.mod
rads_devel_netcdf.o:	rads_devel.mod ../src/rads.mod ../src/rads_misc.mod ../src/rads_netcdf.mod
rads_gen.o:	../src/rads.mod ../src/rads_misc.mod ../src/rads_time.mod
$(PROGS):	rads_devel.o $(LIB)

#-----------------------------------------------------------------------
# Test programs
#-----------------------------------------------------------------------

test:

#-----------------------------------------------------------------------
# How to cleanup
#-----------------------------------------------------------------------

clean spotless:
	$(RM) $(PROGS) *.o *.mod
