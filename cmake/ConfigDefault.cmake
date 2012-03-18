# 
#	$Id: ConfigDefault.cmake 9213 2011-10-06 01:57:06Z fwobbe $
#
#	Copyright (c) 1991-2011 by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis, and F. Wobbe
#	See LICENSE.TXT file for copying and redistribution conditions.
#
#	This program is free software; you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation; version 2 or any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	Contact info: gmt.soest.hawaii.edu
#-------------------------------------------------------------------------------
#
# Useful CMake variables.
#
# There are two configuration files:
#   1) "ConfigDefault.cmake" - is version controlled and used to add new default
#      variables and set defaults for everyone.
#   2) "ConfigUser.cmake" - is not version controlled (currently listed in
#      svn:ignore property) and used to override defaults on a per-user basis.
#
# NOTE: If you want to change CMake behaviour just for yourself then copy
#      "ConfigUserTemplate.cmake" to "ConfigUser.cmake" and then edit
#      "ConfigUser.cmake" (not "ConfigDefault.cmake" or "ConfigUserTemplate.cmake").
#

# The GMT package name.
set (RADS_PACKAGE_NAME "rads4")

# a short description of the gmt project (only a few words).
set (RADS_PACKAGE_DESCRIPTION_SUMMARY "Radar Altimeter Database System 4")

# The GMT package version.
set (RADS_PACKAGE_VERSION_MAJOR "4")
set (RADS_PACKAGE_VERSION_MINOR "0")
set (RADS_PACKAGE_VERSION_PATCH "0")

# The RADS4 package version.
set (RADS_PACKAGE_VERSION "${RADS_PACKAGE_VERSION_MAJOR}.${RADS_PACKAGE_VERSION_MINOR}.${RADS_PACKAGE_VERSION_PATCH}")

# prefer shared libs over static
set (BUILD_SHARED_LIBS true)
set (CMAKE_FIND_STATIC LAST)
