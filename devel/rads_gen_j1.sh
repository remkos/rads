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
# Convert Jason-1 (I)GDR files to RADS
#
# syntax: rads_gen_j1.sh <directories>
#-----------------------------------------------------------------------
. rads_sandbox.sh

rads_open_sandbox j1

date													>  "$log" 2>&1

for tar in "$@"; do
	case "$tar" in
		*.txz) tar -xJf "$tar"; dir=`basename "$tar" .txz` ;;
		*.tgz) tar -xzf "$tar"; dir=`basename "$tar" .tgz` ;;
		*) dir="$tar" ;;
	esac
	ls "$dir"/JA1_???_2P*.nc > "$lst"
	rads_gen_jason	$options < "$lst"					>> "$log" 2>&1
	case "$tar" in
		*.t?z) chmod -R u+w "$dir"; rm -rf "$dir" ;;
	esac
done

# Do the patches to all data

rads_fix_jason    $options --all						>> "$log" 2>&1
rads_add_orbit    $options -Valt_gdre --dir=gdr-e-poe	>> "$log" 2>&1
rads_add_orbit    $options -Valt_gps     -C9-161		>> "$log" 2>&1
rads_add_orbit    $options -Valt_slcci					>> "$log" 2>&1
rads_add_orbit    $options -Valt_std2400 -C1-260		>> "$log" 2>&1
rads_add_ssb      $options --ssb=ssb_tran2012			>> "$log" 2>&1
rads_add_common   $options								>> "$log" 2>&1
rads_add_tide     $options --models=fes14				>> "$log" 2>&1
rads_add_dac      $options --ymd=19910101,20160101 -ue  >> "$log" 2>&1
rads_add_dual     $options -l							>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options -x -x mle3          			>> "$log" 2>&1
rads_add_sla      $options -x -x mle3          			>> "$log" 2>&1

date													>> "$log" 2>&1

rads_close_sandbox
