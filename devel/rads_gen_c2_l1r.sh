#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2021  Remko Scharroo
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
# Convert CryoSat-2 L1R files to RADS
#
# syntax: rads_gen_c2_l1r.sh <directories>
#-----------------------------------------------------------------------

. rads_sandbox.sh

rads_open_sandbox c2

date							>  "$log" 2>&1

for tar in "$@"; do
	case "$tar" in
		*.txz) tar -xJf "$tar"; dir=`basename "$tar" .txz` ;;
		*.tgz) tar -xzf "$tar"; dir=`basename "$tar" .tgz` ;;
		*) dir="$tar" ;;
	esac
	find -L "$dir" -name "CS_*_[B-Z]???.nc" -print | sort -r | sort -u -t/ -k3.20,3.34 > "$lst"
	case "$dir" in
		*/c???) cycle="-C"`basename "$dir" | cut -c2-` ;;
		*)      cycle= ;;
	esac

	rads_gen_c2_l1r $options $cycle < "$lst"	>> "$log" 2>&1

	case "$tar" in
		*.t?z) chmod -R u+w "$dir"; rm -rf "$dir" ;;
	esac
done

case "$dir" in
	*FDM*) orbit_opt="-Valt_gdrf --dir=gdr-f-moe" ;;
	*LRM*) orbit_opt="-Valt_gdrf" ;;
esac

rads_add_ncep    $options -gs               >> "$log" 2>&1
rads_fix_c2      $options --all				>> "$log" 2>&1
rads_add_ssb     $options --all				>> "$log" 2>&1
rads_add_orbit   $options $orbit_opt --equator --loc-7 --rate	>> "$log" 2>&1
rads_add_orbit   $options -Valt_eig6c		>> "$log" 2>&1
rads_add_common  $options 					>> "$log" 2>&1
rads_add_ecmwf   $options --all				>> "$log" 2>&1
rads_add_iono    $options --all				>> "$log" 2>&1
rads_add_mog2d   $options					>> "$log" 2>&1
rads_add_ww3_222 $options --all				>> "$log" 2>&1
rads_add_ww3_314 $options --all -C1-36		>> "$log" 2>&1
rads_add_sla     $options                   >> "$log" 2>&1

date										>> "$log" 2>&1

rads_close_sandbox
