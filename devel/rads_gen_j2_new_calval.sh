#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2016  Remko Scharroo
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
# Convert latest Jason-2 OGDR/IGDR files to RADS, only during Cal/Val phase
#
# The most recently updated data in the OGDR/IGDR files
# will be processed. Four directories will be created, two called
# j2.ogdr0 and j2.igdr0 for the (mostly) original OGDRs, two called
# j2.ogdr and j2.igdr which includes post-processing.
#
# syntax: rads_gen_j2_new_calval.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Check arguments

days=3
types="ogdr igdr"
while getopts "iod:" arg; do
	case $arg in
		d) days=$OPTARG ;;
		o) types=ogdr ;;
		i) types=igdr ;;
	esac
done

# Process only OGDR/IGDR data for the last ($days+1) days (including current)

d0=`date -u -v -${days}d +%Y%m%d 2>&1` || d0=`date -u --date="${days} days ago" +%Y%m%d`

# Run this for ogdr or igdr depending on the arguments on the command line.
# Default is doing both!

for type in ${types}; do

rads_open_sandbox j2.${type}0
lst=$SANDBOX/rads_gen_j2_tmp.lst

date											>  $log 2>&1

omrk=${type}/.bookmark
TZ=UTC touch -t ${d0}0000 $omrk
find ${type}/c??? -name "JA2_*.nc" -a -newer $omrk | sort > $lst
rads_gen_jason --ymd=$d0 $options < $lst			>> $log 2>&1

rads_close_sandbox

# Now process do the same again, and do the post-processing

rads_open_sandbox j2.${type}
find ${type}/c??? -name "JA2_*.nc" -a -newer $omrk | sort > $lst
rads_gen_jason --ymd=$d0 $options < $lst			>> $log 2>&1

# Do the patches to all data

rads_fix_jason   $options --all	--rad=-2.5		>> $log 2>&1
rads_add_ssb     $options --ssb=ssb_tran2012	>> $log 2>&1
rads_add_orbit   $options -Valt_cnes --dir=gdr-e-moe --equator --loc-7 --rate	>> $log 2>&1
rads_add_iono    $options --all					>> $log 2>&1
rads_add_common  $options						>> $log 2>&1
rads_add_dual    $options						>> $log 2>&1
rads_add_dual    $options --ext=mle3			>> $log 2>&1
rads_add_ib      $options						>> $log 2>&1
rads_add_ww3_222 $options --all					>> $log 2>&1
rads_add_sla     $options						>> $log 2>&1
rads_add_sla     $options --ext=mle3			>> $log 2>&1

date											>> $log 2>&1

rads_close_sandbox

done
