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
# Convert latest SWOT files to RADS
#
# The most recently updated data in the directories
# will be processed
#
# syntax: rads_gen_sw_new.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

rads_open_sandbox sw
lst=$SANDBOX/rads_gen_sw_new.lst
imrk=igdr/.bookmark
omrk=ogdr/.bookmark

date																	>  "$log" 2>&1

# Process only OGDR data for the last three days (including current)
d0=`date -u -v -2d +%Y%m%d 2>&1` || d0=`date -u --date="2 days ago" +%Y%m%d`
TZ=UTC touch -t ${d0}0000 $omrk
find -L ogdr/c[0-8]?? -name "SWOT_*.nc" -a -newer $omrk | sort > "$lst"
rads_gen_swot --ymd=$d0 < "$lst"									>> "$log" 2>&1

# Now process all IGDR data that came in during the last four days (including current)
d0=`date -u -v -3d +%Y%m%d 2>&1` || d0=`date -u --date="3 days ago" +%Y%m%d`
TZ=UTC touch -t ${d0}0000 $imrk
find -L igdr/c??? -name "SWOT_*.nc" -a -newer $imrk | sort > "$lst"
rads_gen_swot < "$lst"											>> "$log" 2>&1

# Do the patches to all data

rads_add_common   $options							>> "$log" 2>&1
rads_add_mfwam    $options --wind					>> "$log" 2>&1
rads_add_iono     $options --all					>> "$log" 2>&1
# To support GDR-G with backward compatibility
grep -q _2Pf $lst || rads_add_tide $options --models=fes14		>> "$log" 2>&1

# Redetermine SSHA
rads_add_refframe $options -x -x mle3				>> "$log" 2>&1
rads_add_sla      $options -x -x mle3				>> "$log" 2>&1

date												>> "$log" 2>&1

rads_close_sandbox
