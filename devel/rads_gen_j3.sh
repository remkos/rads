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
# Convert Jason-3 O/I/GDR files to RADS
#
# syntax: rads_gen_j3.sh [ .<type> ] <directories>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Exit when no directory names are provided
[[ $# -eq 0 ]] && exit

type=
case $1 in
	.*) type=$1; shift ;;
esac

rads_open_sandbox j3${type}

date												>  "$log" 2>&1

for tar in "$@"; do
	case "$tar" in
		*cycle[0-9][0-9][0-9]) dir=${tar/cycle/cycle_}; mv "$tar" "$dir" ;;
		*.txz) tar -xJf "$tar"; dir=`basename "$tar" .txz` ;;
		*.tgz) tar -xzf "$tar"; dir=`basename "$tar" .tgz` ;;
		*) dir="$tar" ;;
	esac
	ls "$dir"/JA3_???_2P*.nc > "$lst"
	rads_gen_jason_gdrf	$options < "$lst"			>> "$log" 2>&1
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

rads_fix_jason    $options --all					>> "$log" 2>&1
rads_add_common   $options							>> "$log" 2>&1
rads_add_mfwam    $options -C107-999 --wind			>> "$log" 2>&1
rads_add_orbit    $options -Valt_gps --dir=jplgpspoe -C0-360	>> "$log" 2>&1
rads_add_orbit    $options -Valt_gps --dir=jplgpsmoe -C360-999	>> "$log" 2>&1
rads_add_orbit    $options -Valt_std2400 -C1-513    >> "$log" 2>&1

if grep -q _2Pg $lst ; then
	# For GDR-G we add the FES2014 model
	rads_add_tide $options --models=fes14			>> "$log" 2>&1
else
	# For GDR-F we add MLE3 support
	extra="-x mle3 $extra"
fi

# Redetermine SSHA
rads_add_refframe $options -x $extra				>> "$log" 2>&1
rads_add_sla      $options -x $extra				>> "$log" 2>&1

date												>> "$log" 2>&1

rads_close_sandbox
