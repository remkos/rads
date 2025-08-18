#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2025  Remko Scharroo
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
# Convert Jason-2 GDR-F files to RADS
#
# syntax: rads_gen_j2.sh <directories>
#-----------------------------------------------------------------------
. rads_sandbox.sh

rads_open_sandbox j2

date													>  "$log" 2>&1

for tar in "$@"; do
	case "$tar" in
		*.txz) tar -xJf "$tar"; dir=`basename "$tar" .txz` ;;
		*.tgz) tar -xzf "$tar"; dir=`basename "$tar" .tgz` ;;
		*) dir="$tar" ;;
	esac
	ls "$dir"/JA2_???_2P*.nc > "$lst"
	rads_gen_jason_gdrf	$options < "$lst"				>> "$log" 2>&1
	case "$tar" in
		*.t?z) chmod -R u+w "$dir"; rm -rf "$dir" ;;
	esac
done

# Do the patches to all data

# rads_fix_jason --range only as long as the correction is not in the product
rads_fix_jason    $options --all --range				>> "$log" 2>&1
rads_add_common   $options								>> "$log" 2>&1
rads_add_orbit    $options -Valt_gps     -C1-327		>> "$log" 2>&1
rads_add_orbit    $options -Valt_std2400 -C1-303		>> "$log" 2>&1
rads_add_orbit    $options -Valt_poeg                   >> "$log" 2>&1
rads_add_tide     $options --models=fes14				>> "$log" 2>&1
rads_add_dac      $options --ymd=19910101,20160101 -ue  >> "$log" 2>&1
rads_add_dual     $options -l							>> "$log" 2>&1
rads_add_dual     $options -l --ext=mle3				>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options -x -x mle3 -x adaptive		>> "$log" 2>&1
rads_add_sla      $options -x -x mle3 -x adaptive -Xgdr_g		>> "$log" 2>&1

date													>> "$log" 2>&1

rads_close_sandbox
