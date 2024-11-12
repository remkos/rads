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
# Convert latest SWOT files to RADS, only during Cal/Val phase
#
# The most recently updated data in the files
# will be processed.
#
# syntax: rads_gen_sw_new_calval.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Check arguments

days=3
types="ogdr igdr gdr"
while getopts "iogd:" arg; do
	case $arg in
		d) days=$OPTARG ;;
		o) types=ogdr ;;
		i) types=igdr ;;
		g) types=gdr ;;
	esac
done

# Process only OGDR/IGDR data for the last ($days+1) days (including current)

d0=`date -u -v -${days}d +%Y%m%d 2>&1` || d0=`date -u --date="${days} days ago" +%Y%m%d`

# Run this for ogdr, igdr, or gdr depending on the arguments on the command line.
# Default is doing all three!

for type in ${types}; do

rads_open_sandbox sw.${type}0
lst=$SANDBOX/rads_gen_sw_tmp.lst

date													>  "$log" 2>&1

omrk=${type}/.bookmark
TZ=UTC touch -t ${d0}0000 $omrk
case $type in
	gdr)
		find -L ${type}/cycle_??? -name "SWOT_*.nc" -a -newer $omrk | sort > "$lst"
		if [ -s "$lst" ]; then
			rads_gen_swot $options < "$lst"		>> "$log" 2>&1
		fi
		;;
	*)
		find -L ${type}/c??? -name "SWOT_*.nc" -a -newer $omrk | sort > "$lst"
		rads_gen_swot --ymd=$d0 $options < "$lst"	>> "$log" 2>&1
		;;
esac

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
rads_add_iono     $options --all					>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options -x -x mle3 $extra		>> "$log" 2>&1
rads_add_sla      $options -x -x mle3 $extra		>> "$log" 2>&1

date												>> "$log" 2>&1

rads_close_sandbox

done
