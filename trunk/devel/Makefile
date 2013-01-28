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

include ../config.mk

# Programs to create
PROGS	= rads_gen_reaper

# RADS library and linked libraries
LIB     = ../src/librads.a $(ALTIM)/lib/rssubs.a
LDLIBS  = $(NETCDF_LIB)

#-----------------------------------------------------------------------
# What to create?
#-----------------------------------------------------------------------

all:	$(PROGS)

install:	all
	$(INSTALL_DIR) $(bindir)
	$(INSTALL_PROGRAM) $(PROGS) $(bindir)

man:	$(MAN)

#-----------------------------------------------------------------------
# Object dependencies
#-----------------------------------------------------------------------

$(PROGS):	rads_devel.o $(LIB)
rads_devel.o:	$(LIB)

#-----------------------------------------------------------------------
# How to cleanup
#-----------------------------------------------------------------------
clean spotless:
	$(RM) $(PROGS) *.o *.mod
