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
# Convert latest Sentinel-3A NRT files to RADS
#
# The most recently updated data in the NRT directory
# will be processed.
#
# syntax: rads_gen_3a_calval_new.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

# For the time being: only process NRT data
types=${types:-nrt}
days=${days:-3}

for type in ${types}; do
	d0=`date -u -v -${days}d +%Y%m%d 2>&1` || d0=`date -u --date="${days} days ago" +%Y%m%d`

	mrk=$type/.bookmark
	TZ=UTC touch -t ${d0}0000 $mrk

	dir=3a.${type}0
	rads_open_sandbox $dir
	find $type/c??? -name "*.nc" -a -newer $mrk | sort > $lst
	date >  $log 2>&1
	rads_gen_s3		$options --ymd=$d0 < $lst			>> $log 2>&1
	rads_close_sandbox

# Now process do the same again, and do the post-processing
	dir=3a.${type}1
	rads_open_sandbox $dir
	find $type/c??? -name "*.nc" -a -newer $mrk | sort > $lst
	date >  $log 2>&1
	rads_gen_s3		$options --ymd=$d0 < $lst			>> $log 2>&1

# Make the remaining fixes
	rads_fix_s3     $options --all						>> $log 2>&1
# Recompute SSB
	rads_add_ssb    $options --ssb=ssb_cls				>> $log 2>&1
	rads_add_ssb    $options --ssb=ssb_cls_c			>> $log 2>&1
	rads_add_ssb    $options --ssb=ssb_cls_plrm			>> $log 2>&1
# Recompute dual freq iono and smooth it
	rads_add_dual   $options --recompute				>> $log 2>&1
	rads_add_dual   $options --recompute --ext=plrm		>> $log 2>&1
# Add MOE (and POE) orbit
	rads_add_orbit  $options -Valt_cnes --dir=moe_doris	>> $log 2>&1
	rads_add_orbit  $options -Valt_cnes --dir=poe		>> $log 2>&1
# General geophysical corrections
	rads_add_common $options							>> $log 2>&1
	rads_add_mog2d  $options							>> $log 2>&1
	rads_add_ncep   $options -gdwi						>> $log 2>&1
	rads_add_ecmwf  $options -dwui						>> $log 2>&1
	rads_add_iono   $options -gn						>> $log 2>&1
# Redetermine SSHA
	rads_add_sla    $options							>> $log 2>&1
	rads_add_sla    $options --ext=plrm					>> $log 2>&1

	date

	rads_close_sandbox

done
