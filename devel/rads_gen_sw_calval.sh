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
# Convert SWOT files to RADS
#
# The data from <directory>/<type>/<cycle(s)> will be processed and put
# into RADS.
# Files will be created in
# $RADSDATAROOT/sw.<type>0 - Unadultered files
# $RADSDATAROOT/sw.<type> - RADSified files
#
# syntax: rads_gen_sw_calval.sh <directory>/<types>/<cycles(s)>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Do only latest files when using -d<days>
d0=20000101
case $1 in
	-d*) days=${1:2}; shift
		d0=$(date -u -v -${days}d +%Y%m%d 2>/dev/null || date -u --date="${days} days ago" +%Y%m%d)
		;;
esac

# Exit when no directory names are provided
[[ $# -eq 0 ]] && exit

# Determine type
dir=$(dirname $1)
case $dir in
	*ogdr*) type=ogdr ;;
	*igdr*) type=igdr ;;
	*) type=gdr ;;
esac

# Process "unadultered" files
rads_open_sandbox "sw.${type}0"

date												>  "$log" 2>&1

# Set bookmark according to $d0
mrk=$RADSDATAROOT/.bookmark
TZ=UTC touch -t ${d0}0000 "$mrk"
dirs=
for dir in "$@"; do
	case "$dir" in
		*.txz) tar -C$SANDBOX -xJf "$dir"; dir=$SANDBOX/`basename "$dir" .txz`; chmod -R u+w "$dir" ;;
	esac
	dirs="$dirs $dir"
done
find $dirs -name "SWOT*.nc" -a -newer "$mrk" | sort > "$lst"

# Convert only to RADS, nothing else
rads_gen_swot $options < "$lst"						>> "$log" 2>&1

date												>> "$log" 2>&1

# Now continue with the post-processing
rads_reuse_sandbox "sw.${type}"

# Add adaptive retracker for NTC

case $type in
	gdr) extra="-x adaptive" ;;
	  *) extra= ;;
esac

# Do the patches to all data

rads_add_common   $options							>> "$log" 2>&1
rads_add_mfwam    $options --wind					>> "$log" 2>&1

# If not GDR-F, add the FES2014 model
if grep -q _2Pf $lst || rads_add_tide $options --models=fes14	>> "$log" 2>&1
# If not GDR-G, add MLE3 support
if grep -q _2Pg $lst || extra="-x mle3 $extra"

# Redetermine SSHA
rads_add_refframe $options -x $extra				>> "$log" 2>&1
rads_add_sla      $options -x $extra				>> "$log" 2>&1

date												>> "$log" 2>&1

rads_close_sandbox
