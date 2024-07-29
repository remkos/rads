#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2024  Remko Scharroo
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
#
# Add/overwrite numerous fields in the RADS data base related to GDR-G standards (#204):
# - mss_dtu21 (#201, #174)
# - mss_hybrid23 (#201)
# - fes22,lptide (#200)
# and some other needed patches as per the following issues:
# - srtm15plus (#203)
# - oisst v2.1 (#187)
# - ptide (#176)
# - era5 (#142)
#
rads_add_tide     "$@" --models=fes22,lptide,ptide
rads_add_grid     "$@" -Vmss_dtu21,mss_hybrid23,topo_srtm15plus
rads_add_sst      "$@" --all
rads_add_era5     "$@" --all --new
