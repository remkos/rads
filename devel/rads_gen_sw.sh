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
# Convert SWOT files to RADS
#
# syntax: rads_gen_sw.sh [ .<type> ] <directories>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Exit when no directory names are provided
[[ $# -eq 0 ]] && exit

type=
case $1 in
	.*) type=$1; shift ;;
esac

rads_open_sandbox sw${type}
lst=$SANDBOX/rads_gen_sw.lst

date												>  "$log" 2>&1

for tar in "$@"; do
	case "$tar" in
		*cycle[0-9][0-9][0-9]) dir=${tar/cycle/cycle_}; mv "$tar" "$dir" ;;
		*.txz) tar -xJf "$tar"; dir=`basename "$tar" .txz` ;;
		*.tgz) tar -xzf "$tar"; dir=`basename "$tar" .tgz` ;;
		*) dir="$tar" ;;
	esac
	ls "$dir"/SWOT_???_2Pf*.nc > "$lst"
	rads_gen_swot	$options < "$lst"			>> "$log" 2>&1
	case "$tar" in
		*.t?z) chmod -R u+w "$dir"; rm -rf "$dir" ;;
	esac
done

# Add adaptive retracker for NTC

case $type in
	ogdr|igdr) extra= ;;
	        *) extra="-x adaptive" ;;
esac

# Do the patches to all data

rads_add_common   $options							>> "$log" 2>&1
rads_add_mfwam    $options --all					>> "$log" 2>&1
rads_add_iono     $options --all					>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options -x -x mle3 $extra		>> "$log" 2>&1
rads_add_sla      $options -x -x mle3 $extra		>> "$log" 2>&1

date												>> "$log" 2>&1

rads_close_sandbox
