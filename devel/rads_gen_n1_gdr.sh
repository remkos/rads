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
# Convert Envisat GDR v3 files to RADS
#
# syntax: rads_gen_n1_gdr.sh <directories>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Exit when no directory names are provided
[[ $# -eq 0 ]] && exit

rads_open_sandbox n1.gdr

find "$@" -name "*.nc"| sort > "$lst"
date													>  "$log" 2>&1
rads_gen_n1_gdr  $options < "$lst"						>> "$log" 2>&1

# Add the FDR4ALT data

find $RADSROOT/ext/FDR4ALT/TDP_??/en1/c??? -name "*.nc" | sort -t_ -k7,8 > "$lst"
date
rads_add_f4a     $options < "$lst"						>> "$log" 2>&1

# Do the patches to all data

rads_add_iono    $options --all							>> "$log" 2>&1
rads_add_common  $options								>> "$log" 2>&1
rads_add_era5    $options --all							>> "$log" 2>&1
rads_fix_n1      $options --all							>> "$log" 2>&1

# Redetermine SSHA
rads_add_refframe $options -x -x adaptive				>> "$log" 2>&1
rads_add_sla      $options -x -x adaptive				>> "$log" 2>&1

date													>> "$log" 2>&1

rads_close_sandbox
