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
# Convert latest Jason-3 OGDR/IGDR/GDR files to RADS, only during Cal/Val phase
#
# The most recently updated data in the OGDR/IGDR/GDR files
# will be processed. Four directories will be created, three called
# j3.ogdr0, j3.igdr0, and j3.gdr0 for the (mostly) original O/I/GDRs, three called
# j3.ogdr, j3.igdr, and j3.gdr which include post-processing.
#
# syntax: rads_gen_j3_new_calval.sh
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

rads_open_sandbox j3.${type}0
lst=$SANDBOX/rads_gen_j3_tmp.lst

date												>  "$log" 2>&1

omrk=${type}/.bookmark
TZ=UTC touch -t ${d0}0000 $omrk
case $type in
	gdr)
		find -L ${type}/cycle_??? -name "JA3_*.nc" -a -newer $omrk | sort > "$lst"
		if [ -s "$lst" ]; then
			rads_gen_jason_gdrf $options < "$lst"			>> "$log" 2>&1
		fi
		;;
	*)
		find -L ${type}/c??? -name "JA3_*.nc" -a -newer $omrk | sort > "$lst"
		rads_gen_jason_gdrf --ymd=$d0 $options < "$lst"	>> "$log" 2>&1
		;;
esac

rads_close_sandbox

# Now process do the same again, and do the post-processing

rads_open_sandbox j3.${type}
case $type in
	gdr)
		find -L ${type}/cycle_??? -name "JA3_*.nc" -a -newer $omrk | sort > "$lst"
		if [ -s "$lst" ]; then
			rads_gen_jason_gdrf $options < "$lst"			>> "$log" 2>&1
		fi
		;;
	*)
		find -L ${type}/c??? -name "JA3_*.nc" -a -newer $omrk | sort > "$lst"
		rads_gen_jason_gdrf --ymd=$d0 $options < "$lst"	>> "$log" 2>&1
		;;
esac

# Do the patches to all data

rads_fix_jason    $options --all					>> "$log" 2>&1
rads_add_common   $options							>> "$log" 2>&1
rads_add_iono     $options --all					>> "$log" 2>&1
rads_add_sla      $options							>> "$log" 2>&1
rads_add_refframe $options --ext=mle3				>> "$log" 2>&1
rads_add_sla      $options --ext=mle3				>> "$log" 2>&1
case $type in
	gdr)
		rads_add_refframe  $options --ext=adaptive	>> "$log" 2>&1
		rads_add_sla       $options --ext=adaptive	>> "$log" 2>&1
		;;
esac

date												>> "$log" 2>&1

rads_close_sandbox

done
