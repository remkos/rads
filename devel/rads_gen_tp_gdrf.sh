#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2023  Remko Scharroo
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
# Convert TOPEX/Poseidon Retracked GDR-F files to RADS
#
# syntax: rads_gen_tp_gdrf.sh <directories>
#-----------------------------------------------------------------------
sat=$1
shift

. rads_sandbox.sh

rads_open_sandbox ${sat}.gdrf
lst=$SANDBOX/rads_gen_${sat}_gdrf.lst

date													>  "$log" 2>&1

for tar in "$@"; do
	case "$tar" in
		*.txz) tar -xJf "$tar"; dir=`basename "$tar" .txz` ;;
		*.tgz) tar -xzf "$tar"; dir=`basename "$tar" .tgz` ;;
		*) dir="$tar" ;;
	esac
	ls "$dir"/TP_GPN_*.nc > "$lst"
	rads_gen_tp_gdrf $options < "$lst"					>> "$log" 2>&1
	case "$tar" in
		*.t?z) chmod -R u+w "$dir"; rm -rf "$dir" ;;
	esac
done

# Do the patches to all data

# GIM iono is only available from 1994-01-01 onward
rads_add_iono     $options -C1-47 --nic09				>> "$log" 2>&1
rads_add_iono     $options -C48-481 --all				>> "$log" 2>&1

rads_add_common   $options								>> "$log" 2>&1
case $sat in
tx) rads_add_dual $options -r							>> "$log" 2>&1
    rads_add_dual $options -r -x mle3					>> "$log" 2>&1
    extra="-x mle3" ;;
esac

# Redetermine SSHA and add ERA5 (temporarily suppressed)
#rads_add_refframe $options -x $extra					>> "$log" 2>&1
#rads_add_sla      $options -x $extra					>> "$log" 2>&1
#rads_add_era5     $options --all						>> "$log" 2>&1

date													>> "$log" 2>&1

rads_close_sandbox
