#-----------------------------------------------------------------------
# $Id$
#
# Copyright (C) 2011  Remko Scharroo (Altimetrics LLC)
# See LICENSE.TXT file for copying and redistribution conditions.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#-----------------------------------------------------------------------
#
# To build out-of-source do (example):
#
#	mkdir build
#	cd build
#	cmake -DCMAKE_INSTALL_PREFIX=gmt5 -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
#
# CMAKE_BUILD_TYPE can be: empty, Debug, Release, RelWithDebInfo or MinSizeRel
#
# cmake creates a new file cmake/ConfigUser.cmake if it does not already
# exist. You can configure additional options there.
#

# Specify minimum required cmake version
cmake_minimum_required(VERSION 2.8)

# Define project name
project (RADS4 Fortran)

# Where to find our CMake modules (this variable is visible in subdirectories).
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Add subdirectories
#add_subdirectory (doc)
add_subdirectory (src)

# Configuration done
message (
	"*  Installing RADS in   : ${CMAKE_INSTALL_PREFIX}\n"
	"*  RADSDATADIR          : ${RADSDATADIR}\n"
	"*  NetCDF library       : ${NETCDF_LIBRARY}\n"
	"*  NetCDF include dir   : ${NETCDF_INCLUDE_DIR}\n"
)
