#-----------------------------------------------------------------------
# config.mk.  Generated from config.mk.in by configure.
#
# Copyright (c) 2011-2022  Remko Scharroo
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

#-----------------------------------------------------------------------
# Package information
#-----------------------------------------------------------------------
PACKAGE_TARNAME = rads
PACKAGE_VERSION = 4.4.0

#-----------------------------------------------------------------------
# Object locations: Library, Manuals, Binaries
#-----------------------------------------------------------------------
prefix = /opt/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
datadir = /Users/davidtrossman/Documents/Code/rads/data
bindir = ${exec_prefix}/bin
docdir = ${datarootdir}/doc/${PACKAGE_TARNAME}
includedir = ${prefix}/include
libdir = ${exec_prefix}/lib
RADS_LIB = -L$(libdir) -lrads
RADS_INC = -I$(includedir)
INSTALL = /opt/local/bin/ginstall -c
INSTALL_PROGRAM = ${INSTALL} -C
INSTALL_SCRIPT = ${INSTALL} -C
INSTALL_DATA = ${INSTALL} -m 644 -C
INSTALL_DIR = ${INSTALL} -d

#-----------------------------------------------------------------------
# Are we little endian?
#-----------------------------------------------------------------------
LITTLE_ENDIAN = true

#-----------------------------------------------------------------------
# Environments for developers only
#-----------------------------------------------------------------------
ALTIM = /Users/davidtrossman/Documents/Code/rads/altim
DEVEL_DIRS = devel
DEVEL_OBJS = ${ALTIM}/src/rssubs/tpj_subs.o ${ALTIM}/src/rssubs/solar_subs.o
LAPACK_LIB = 

#-----------------------------------------------------------------------
# Include and load option for NetCDF
#-----------------------------------------------------------------------
NETCDF_INC = -I/opt/local/include/
NETCDF_LIB = -L/opt/local/lib/ -lnetcdff -lnetcdf

#-----------------------------------------------------------------------
# Fortran compiler settings
#-----------------------------------------------------------------------
FC = /opt/local/bin/gfortran-mp-12
FFLAGS = -O2 $(INCLUDES)
MOD_FLAG = -I
DEVEL_FLAGS =  -I../devel -I${ALTIM}/include

#-----------------------------------------------------------------------
# Additional tools
#-----------------------------------------------------------------------
LN_S = ln -s -f
RANLIB = ranlib

#-----------------------------------------------------------------------
# Special rules
#-----------------------------------------------------------------------
# .f90 rules
%.o:	%.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $<
%:	%.f90
	$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

# Cancel default .mod rules
%:	%.mod
%.o:	%.mod

# Create new .mod rule
# Do not add %.f90 dependency because %.mod retains timestamp when not altered
%.mod:
	$(COMPILE.f) -o $*.o $*.f90

# In case of Git only
GIT_INDEX	= @GIT_INDEX@
