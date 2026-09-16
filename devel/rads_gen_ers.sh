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

# Make a temporary list file
tmplst=$(mktemp ${TMPDIR:-/tmp}/rads.XXXXXX)
find "$@" -name "*.nc"| sort > "$tmplst"

# If list is empty, exit
[[ ! -s $tmplst ]] && rm -f "$tmplst" && exit

# Determine which satellite
sat=e1
grep -q E2_REAP_ERS "$tmplst" && sat=e2

rads_open_sandbox $sat.f4a

mv -f "$tmplst" "$lst"

date													>  "$log" 2>&1
rads_gen_reaper   $options < "$lst"						>> "$log" 2>&1

# Add the FDR4ALT data. Each type will need to be added one at a time

date													>> "$log" 2>&1
for tdp in TDP_OC TDP_WA TDPATM ; do
	find $RADSROOT/ext/FDR4ALT/$tdp/er${sat:1:1}/c??? -name "*.nc" | sort > "$lst"
	rads_add_f4a_ers  $options < "$lst"					>> "$log" 2>&1
done

# Do the patches to all data

date													>> "$log" 2>&1
rads_add_tide     $options --models=fes14				>> "$log" 2>&1
rads_add_common   $options								>> "$log" 2>&1

# Redetermine SSHA
rads_add_refframe $options								>> "$log" 2>&1
rads_add_sla      $options								>> "$log" 2>&1

date													>> "$log" 2>&1

rads_close_sandbox
