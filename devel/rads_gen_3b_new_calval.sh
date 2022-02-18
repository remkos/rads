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
# Convert latest Sentinel-3B NRT and STC files to RADS
#
# The most recently updated data in the NRT and STC directories
# will be processed.
#
# syntax: rads_gen_3b_calval_new.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Process NRT/STC/NTC data
days=${days:-3}
types="${types:-nrt stc}"
while getopts "nsd:" arg; do
	case $arg in
		d) days=$OPTARG ;;
		n) types=nrt ;;
		s) types=stc ;;
	esac
done

d0=`date -u -v -${days}d +%Y%m%d 2>&1` || d0=`date -u --date="${days} days ago" +%Y%m%d`

for type in ${types}; do
	mrk=$type/.bookmark
	TZ=UTC touch -t ${d0}0000 "$mrk"

	dir=3b.${type}0
	rads_open_sandbox "$dir"
	find $type/c??? -name "*.nc" -a -newer "$mrk" | sort > "$lst"
	date >  "$log" 2>&1
	rads_gen_s3		$options --min-rec=6 --ymd=$d0 < "$lst"	>> "$log" 2>&1
	rads_close_sandbox

# Now process do the same again, and do the post-processing
	dir=3b.${type}1
	rads_open_sandbox "$dir"
	find $type/c??? -name "*.nc" -a -newer "$mrk" | sort > "$lst"
	date >  "$log" 2>&1
	rads_gen_s3		$options --min-rec=6 --ymd=$d0 < "$lst"	>> "$log" 2>&1
	rads_fix_s3		$options --all							>> "$log" 2>&1

# Add MOE orbit (for NRT and STC only)
	case $type in
		nr*|st*) rads_add_orbit  $options -Valt_cnes --dir=moe_doris	>> "$log" 2>&1
	esac

# General geophysical corrections
    rads_add_grid     $options -Vangle_coast                >> "$log" 2>&1
	rads_add_common   $options								>> "$log" 2>&1
	rads_add_mfwam    $options --all						>> "$log" 2>&1
	rads_add_iono     $options --all						>> "$log" 2>&1
# Redetermine SSHA
	rads_add_refframe $options -x -x plrm					>> "$log" 2>&1
	rads_add_sla      $options -x -x plrm					>> "$log" 2>&1

	date													>> "$log" 2>&1

	rads_close_sandbox

done
