#!/bin/bash
#-----------------------------------------------------------------------
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
#
# Convert Sentinel-3A files to RADS
#
# The data from <directory>/<type>/<cycle(s)> will be processed and put
# into RADS.
# Files will be created in
# $RADSDATAROOT/3a.<type>0 - Unadultered files
# $RADSDATAROOT/3a.<type>1 - RADSified files
#
# syntax: rads_gen_3a_calval.sh <directory>/<type>/<cycles(s)>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Exit when no directory names are provided
[[ $# -eq 0 ]] && exit

# Determine type
type=$(dirname $1)
type=$(basename $type)

# Process "unadultered" files
rads_open_sandbox "3a.${type}0"

date													>  "$log" 2>&1

find "$@" -name "S3A_*.nc" -o -name "standard_measurement.nc" | sort	> "$lst"
rads_gen_s3 	$options --min-rec=6 < "$lst"			>> "$log" 2>&1

# Now continue with the post-processing
rads_reuse_sandbox "3a.${type}1"

date													>> "$log" 2>&1

# General geophysical corrections
rads_add_grid     $options -Vangle_coast                >> "$log" 2>&1
rads_add_common   $options								>> "$log" 2>&1
rads_add_mfwam    $options -C40-199 --all				>> "$log" 2>&1
rads_add_iono     $options --all						>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options -x -x plrm					>> "$log" 2>&1
rads_add_sla      $options -x -x plrm					>> "$log" 2>&1

date													>> "$log" 2>&1

rads_close_sandbox
