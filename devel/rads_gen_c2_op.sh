#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2025  Remko Scharroo and Eric Leuliette
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
# Convert latest CryoSat-2 NOP, IOP, and GOP files to RADS
#
# The most recently updated data in the nop, iop, and gop directories
# will be processed.
#
# syntax: rads_gen_c2_op.sh <directories>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Exit when no directory names are provided
[[ $# -eq 0 ]] && exit

rads_open_sandbox c2

find "$@" -name "*.nc"| sort > "$lst"
date >  "$log" 2>&1
rads_gen_c2_op		$options < "$lst"	>> "$log" 2>&1

# General geophysical corrections
rads_add_common   $options				>> "$log" 2>&1
rads_add_iono     $options --all		>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options				>> "$log" 2>&1
rads_add_sla      $options				>> "$log" 2>&1

date									>> "$log" 2>&1

rads_close_sandbox
