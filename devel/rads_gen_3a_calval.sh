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
# Convert Sentinel-3A files to RADS
#
# The data from <directory>/<type>/<cycle(s)> will be processed and put
# into RADS.
# Files will be created in
# $RADSDATAROOT/3a.<type>0 - Unadultered files
# $RADSDATAROOT/3a.<type>1 - RADSified files
#
# syntax: rads_gen_3a_calval.sh <directory>/<types>/<cycles(s)>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Determine type
type=$(dirname $1)
type=$(basename $type)

# Process "unadultered" files
org=3a.${type}0
rads_open_sandbox $org
find $* -name standard_measurement.nc | sort	> $lst
date >  $log 2>&1
rads_gen_s3 -S${org} --ymd=$d0 < $lst				>> $log 2>&1
rads_close_sandbox

# Now process do the same again, and do the post-processing
fix=3a.${types}1
rads_open_sandbox $fix
find $* -name standard_measurement.nc | sort	> $lst
date >  $log 2>&1
rads_gen_s3 -S${fix} --ymd=$d0 < $lst				>> $log 2>&1

# Make the remaining fixes
rads_fix_s3     -S${fix} --all						>> $log 2>&1
# Recompute SSB
rads_add_ssb    -S${fix} --ssb=ssb_cls				>> $log 2>&1
rads_add_ssb    -S${fix} --ssb=ssb_cls_c			>> $log 2>&1
rads_add_ssb    -S${fix} --ssb=ssb_cls_plrm			>> $log 2>&1
# Recompute dual freq iono and smooth it
rads_add_dual   -S${fix} --recompute				>> $log 2>&1
rads_add_dual   -S${fix} --recompute --ext=plrm		>> $log 2>&1
# Add MOE (and POE) orbit
rads_add_orbit  -S${fix} -Valt_cnes --dir=moe_doris	>> $log 2>&1
rads_add_orbit  -S${fix} -Valt_cnes --dir=poe		>> $log 2>&1
# General geophysical corrections
rads_add_common -S${fix}							>> $log 2>&1
rads_add_mog2d  -S${fix}							>> $log 2>&1
rads_add_ncep   -S${fix} -gdwi						>> $log 2>&1
rads_add_ecmwf  -S${fix} -dwui						>> $log 2>&1
rads_add_iono   -S${fix} -gn						>> $log 2>&1
# Redetermine SSHA
rads_add_sla    -S${fix}							>> $log 2>&1
rads_add_sla    -S${fix} --ext=plrm					>> $log 2>&1

date												>> $log 2>&1
rads_close_sandbox
