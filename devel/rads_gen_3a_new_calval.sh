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
days=${days:-2}
types="${types:-nrt}"

# Process only NRT data for the last three days (including current)
d0=`date -u -v -${days}d +%Y%m%d 2>&1` || d0=`date -u --date="${days} days ago" +%Y%m%d`

for type in ${types}; do

	mrk=$type/.bookmark
	TZ=UTC touch -t ${d0}0000 $mrk

	org=3a.${type}0
	rads_open_sandbox $org
	lst=$SANDBOX/rads_gen_3a_new_calval.lst
	find $type/c??? -name "*.nc" -a -newer $mrk | sort > $lst
	date >  $log 2>&1
	rads_gen_s3 -S${org} --ymd=$d0 < $lst >> $log 2>&1
	rads_close_sandbox

# Now process do the same again, and do the post-processing
	fix=3a.${types}1
	rads_open_sandbox $fix
	lst=$SANDBOX/rads_gen_3a_new_calval.lst
	find $type/c??? -name "*.nc" -a -newer $mrk | sort > $lst
	date >  $log 2>&1
	rads_gen_s3 -S${fix} --ymd=$d0 < $lst >> $log 2>&1

# Make the remaining fixes
	rads_fix_s3     -S${fix} --all >> $log 2>&1
# Recompute SSB
	rads_add_ssb    -S${fix} --ssb=ssb_cls >> $log 2>&1
	rads_add_ssb    -S${fix} --ssb=ssb_cls_c >> $log 2>&1
	rads_add_ssb    -S${fix} --ssb=ssb_cls_plrm >> $log 2>&1
# Recompute dual freq iono and smooth it
	rads_add_dual   -S${fix} --recompute >> $log 2>&1
	rads_add_dual   -S${fix} --recompute --ext=plrm >> $log 2>&1
# Add MOE (and POE) orbit
	rads_add_orbit  -S${fix} -Valt_cnes --dir=moe_doris >> $log 2>&1
	rads_add_orbit  -S${fix} -Valt_cnes --dir=poe >> $log 2>&1
# General geophysical corrections
	rads_add_common -S${fix} >> $log 2>&1
	rads_add_mog2d  -S${fix} >> $log 2>&1
	rads_add_ncep   -S${fix} -gdwi >> $log 2>&1
	rads_add_iono   -S${fix} -gn >> $log 2>&1
# Redetermine SSHA
	rads_add_sla    -S${fix} >> $log 2>&1
	rads_add_sla    -S${fix} --ext=plrm >> $log 2>&1

	date

	rads_close_sandbox

done
