#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2018  Remko Scharroo
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
dir=3a.${type}0
rads_open_sandbox $dir

date												>  $log 2>&1

find $* -name "*.nc" | sort		> $lst
rads_gen_s3 	$options --min-rec=6 < $lst			>> $log 2>&1
rads_close_sandbox

# Now process do the same again, and do the post-processing
dir=3a.${type}1
rads_open_sandbox $dir

date												>  $log 2>&1

find $* -name "*.nc" | sort		> $lst
rads_gen_s3 	$options --min-rec=6 < $lst			>> $log 2>&1

# Make the remaining fixes
rads_fix_s3     $options --all						>> $log 2>&1
# Recompute SSB
rads_add_ssb    $options --ssb=ssb_cls				>> $log 2>&1
rads_add_ssb    $options --ssb=ssb_cls_c			>> $log 2>&1
rads_add_ssb    $options --ssb=ssb_cls_plrm			>> $log 2>&1
# Recompute dual freq iono and smooth it
rads_add_dual   $options --recompute				>> $log 2>&1
rads_add_dual   $options --recompute --ext=plrm		>> $log 2>&1
# General geophysical corrections
rads_add_common $options							>> $log 2>&1
rads_add_mog2d  $options							>> $log 2>&1
rads_add_ncep   $options -gdwi						>> $log 2>&1
rads_add_ecmwf  $options -dwui						>> $log 2>&1
rads_add_iono   $options --all						>> $log 2>&1
# Redetermine SSHA
rads_add_sla    $options							>> $log 2>&1
rads_add_sla    $options --ext=plrm					>> $log 2>&1

date												>> $log 2>&1

rads_close_sandbox
