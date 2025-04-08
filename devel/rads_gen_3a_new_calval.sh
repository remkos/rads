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
# Convert latest Sentinel-3 NRT and STC files to RADS
#
# The most recently updated data in the NRT and STC directories
# will be processed.
#
# syntax: rads_gen_3?_calval_new.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Which satellite?
case $0 in
	*rads_gen_3a*) sat=3a ;;
	*rads_gen_3b*) sat=3b ;;
	*) echo "$0: unknown script" ; exit ;;
esac

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

d0=$(date -u -v -${days}d +%Y%m%d 2>/dev/null || date -u --date="${days} days ago" +%Y%m%d)

for type in ${types}; do
	rads_open_sandbox "${sat}.${type}0"

	date													>  "$log" 2>&1
	mrk=$RADSDATAROOT/.bookmark
	TZ=UTC touch -t ${d0}0000 "$mrk"
	find $type/c??? -name "*.nc" -a -newer "$mrk" | sort > "$lst"
	rads_gen_s3		$options --min-rec=6 --ymd=$d0 < "$lst"	>> "$log" 2>&1

# Now continue with the post-processing
	rads_reuse_sandbox "${sat}.${type}1"

	date													>> "$log" 2>&1

# Add MOE orbit (for NRT and STC only) and CPOD POE (for NTC/REP only)
	case $type in
		nr*|st*) rads_add_orbit  $options -Valt_cnes --dir=moe_doris	>> "$log" 2>&1 ;;
		*)	     rads_add_orbit  $options -Valt_cpod					>> "$log" 2>&1 ;;
	esac

# General geophysical corrections
	rads_add_common   $options								>> "$log" 2>&1
	rads_add_mfwam    $options --all --new					>> "$log" 2>&1
	rads_add_iono     $options --all						>> "$log" 2>&1
# Redetermine SSHA
	rads_add_refframe $options -x -x plrm					>> "$log" 2>&1
	rads_add_sla      $options -x -x plrm					>> "$log" 2>&1

	date													>> "$log" 2>&1

	rads_close_sandbox

done
