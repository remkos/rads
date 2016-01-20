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
# Convert latest Jason-3 OGDR files to RADS, only during Cal/Val phase
#
# The most recently updated data in the OGDR directories
# will be processed. Two directories will be created, one called
# j3.ogdr_orig for the (mostly) original OGDRs, one called j3.ogdr
# which includes post-processing.
#
# syntax: rads_gen_j3_new.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

rads_open_sandbox j3.ogdr0
lst=$SANDBOX/rads_gen_j3_tmp.lst
omrk=ogdr/.bookmark

date								>  $log 2>&1

# Process only OGDR data for the last three days (including current)

d0=`date -u -v -2d +%Y%m%d 2>&1` || d0=`date -u --date="2 days ago" +%Y%m%d`
TZ=UTC touch -t ${d0}0000 $omrk
find ogdr/c??? -name "JA3_*.nc" -a -newer $omrk | sort > $lst
rads_gen_j3 --ymd=$d0 $options < $lst		>> $log 2>&1

rads_close_sandbox

# Now process do the same again, and do the post-processing

rads_open_sandbox j3.ogdr
find ogdr/c??? -name "JA3_*.nc" -a -newer $omrk | sort > $lst
rads_gen_j3 --ymd=$d0 $options < $lst		>> $log 2>&1

# Do the patches to all data

rads_fix_j3      $options --all		>> $log 2>&1
rads_add_ssb     $options --ssb=ssb_tran2012	>> $log 2>&1
rads_add_orbit   $options -Valt_gdre --dir=gdr-e-moe --equator --loc-7 --rate	>> $log 2>&1
rads_add_iono    $options --all		>> $log 2>&1
rads_add_common  $options			>> $log 2>&1
rads_add_dual    $options			>> $log 2>&1
rads_add_dual    $options --mle=3	>> $log 2>&1
rads_add_ib      $options			>> $log 2>&1
rads_add_ww3_222 $options --all		>> $log 2>&1
rads_add_sla     $options           >> $log 2>&1

date								>> $log 2>&1

rads_close_sandbox
