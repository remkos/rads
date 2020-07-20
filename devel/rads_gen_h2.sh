#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2020  Remko Scharroo
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
# Convert HY-2A IGDR files to RADS
#
# syntax: rads_gen_h2.sh <directories>
#-----------------------------------------------------------------------
. rads_sandbox.sh

rads_open_sandbox 2a
lst=$SANDBOX/rads_gen_h2.lst

date											>  "$log" 2>&1

for tar in "$@"; do
	case "$tar" in
		*.txz) tar -xJf "$tar"; dir=`basename "$tar" .txz` ;;
		*.tgz) tar -xzf "$tar"; dir=`basename "$tar" .tgz` ;;
		*) dir="$tar" ;;
	esac
	ls "$dir"/H2A_RA1_IDR_2PT_*.nc > "$lst"
	rads_gen_h2	$options < "$lst"				>> "$log" 2>&1
	case "$tar" in
		*.t?z) chmod -R u+w "$dir"; rm -rf "$dir" ;;
	esac
done

# Do the patches to all data

#rads_fix_h2      $options --all					>> "$log" 2>&1
#rads_add_ssb     $options --ssb=ssb_tran2012	>> "$log" 2>&1
rads_add_iono    $options --all					>> "$log" 2>&1
rads_add_common  $options						>> "$log" 2>&1
rads_add_dual    $options						>> "$log" 2>&1
rads_add_ib      $options						>> "$log" 2>&1
rads_add_ecmwf   $options --all --wind			>> "$log" 2>&1
#rads_add_orbit   $options -Valt_gdre    -C0-253	>> "$log" 2>&1
rads_add_mog2d   $options						>> "$log" 2>&1
rads_add_ww3_222 $options --all					>> "$log" 2>&1
rads_add_sla     $options           			>> "$log" 2>&1

date											>> "$log" 2>&1

rads_close_sandbox
