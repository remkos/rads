#-----------------------------------------------------------------------
# Copyright (c) 2011-2019  Remko Scharroo
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

# Include file
CONFIG  = config.f90

# Objects
OBJS1   = rads_time.o rads_misc.o rads_grid.o rads_geo.o rads_math.o rads_netcdf.o xmlparse.o
OBJS    = $(OBJS1) rads.o rads3.o

# .mod files created by the objects
MODS1   = $(OBJS1:.o=.mod)
MODS    = $(OBJS:.o=.mod)

# Includes
INCLUDES= $(DEVEL_FLAGS) -Iinclude $(MOD_FLAG). $(NETCDF_INC)

# RADS library and linked libraries
LIB     = librads.a
LDLIBS  = $(NETCDF_LIB)

# Programs created by objects
PROGS   = rads2asc rads2adr rads2nc rads2grd radsstat radscolin radsxogen radsxoconv \
	radsxolist rads_test radsreconfig rads_sample rads3test radsvar \
	radspassesindex

# Programs created by linking to other programs
PROGS_L1= rads2xadr rads2xgf
PROGS_L2= radsxostat
PROGS_L4= rads2asc4 radscolin4 radsstat4 rads2grd4 radspassesindex4

# Shell scripts
SCRIPTS = rads-config

#-----------------------------------------------------------------------
# What to create?
#-----------------------------------------------------------------------

all:	$(PROGS) $(PROGS_L1) $(PROGS_L2) $(PROGS_L4)

install:	all
	$(INSTALL_DIR) $(libdir) $(includedir) $(bindir)
	$(INSTALL_DATA) $(LIB) $(libdir)
	$(INSTALL_DATA) $(MODS) $(includedir)
	$(INSTALL_PROGRAM) $(PROGS) $(bindir)
	$(INSTALL_SCRIPT) $(SCRIPTS) $(bindir)
	for p in $(PROGS_L1); do ( cd $(bindir) ; $(LN_S) rads2adr $$p ) ; done
	for p in $(PROGS_L2); do ( cd $(bindir) ; $(LN_S) radsxolist $$p ) ; done
	for p in $(PROGS_L4); do ( cd $(bindir) ; $(LN_S) $${p%4} $$p ) ; done

#-----------------------------------------------------------------------
# Object dependencies
#-----------------------------------------------------------------------

rads_math.o:	rads_misc.mod rads_time.mod
rads_netcdf.o:	rads_misc.mod rads_netcdf_get_var.f90 include/nf90_message.f90 ../devel/nf90_message.f90
rads.o:			$(MODS1) $(CONFIG) include/rads_tpj.f90 ../devel/rads_tpj.f90
rads3.o:		rads.mod
$(PROGS):		$(LIB)
radsxogen:		spline.o

#-----------------------------------------------------------------------
# Create configuration files
#-----------------------------------------------------------------------

config.f90:	../config.mk
	printf "character(len=rads_naml), parameter :: radsdataroot = '%s'\n" $(datadir) > $@
	printf "logical, parameter :: little_endian = .%s.\n" $(LITTLE_ENDIAN) >> $@
	printf "character(len=rads_varl), parameter :: rads_version_id = '%s'\n" $(PACKAGE_VERSION) >> $@

rads-config:	rads-config.in
	@echo "Error: Need to run configure in main directory first"
	@exit 1

#-----------------------------------------------------------------------
# How to create the library
#-----------------------------------------------------------------------

$(LIB):	$(OBJS) $(DEVEL_OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS) $(DEVEL_OBJS)
	$(RANLIB) $@

#-----------------------------------------------------------------------
# Test programs
#-----------------------------------------------------------------------

TEST=rads_test -vvv -Se1/a --ymd=19910924,19910925 -R-90/120/-50/30 -Vtime,lat,lon,sla,12_,tide_ocean_,flags_,pass,drange | egrep -v "Created|CPU" | sed "s@$(datadir)@@g"
test:	rads_test
	$(TEST) | diff - rads_test.out
rads_test.out:	rads_test
	$(TEST) > $@

#-----------------------------------------------------------------------
# Make links to executables
#-----------------------------------------------------------------------

$(PROGS_L1):	rads2adr
	$(LN_S) rads2adr $@

$(PROGS_L2):	radsxolist
	$(LN_S) radsxolist $@

$(PROGS_L4):	%4:	%
	$(LN_S) $* $@

#-----------------------------------------------------------------------
# How to cleanup
#-----------------------------------------------------------------------

clean:
	$(RM) $(LIB) $(PROGS) *.o *.mod
spotless:	clean
	$(RM) $(CONFIG)
