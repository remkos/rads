#!/bin/bash
#-----------------------------------------------------------------------
# $Id$
#
# Copyright (c) 2011-2015  Remko Scharroo
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
# Convert SARAL O/I/GDR files to RADS
#
# syntax: rads_gen_saral.sh <directories>
#-----------------------------------------------------------------------
. rads_sandbox.sh

rads_open_sandbox sa
lst=$SANDBOX/rads_gen_saral.lst

date													>  $log 2>&1

for tar in $*; do
	case $tar in
		*.txz) tar -xJf $tar; dir=`basename $tar .txz` ;;
		*.tgz) tar -xzf $tar; dir=`basename $tar .tgz` ;;
		*) dir=$tar ;;
	esac
	ls $dir/SRL_*.nc > $lst
	rads_gen_saral $options < $lst						>> $log 2>&1
	case $tar in
		*.t?z) chmod -R u+w $dir; rm -rf $dir ;;
	esac
done

# Do the patches to all data

rads_add_iono    $options --all							>> $log 2>&1
rads_add_common  $options								>> $log 2>&1
rads_add_ncep    $options -gdwu --sig0-saral			>> $log 2>&1
rads_fix_sa      $options								>> $log 2>&1
rads_add_ssb     $options --all							>> $log 2>&1
rads_add_mog2d   $options								>> $log 2>&1
rads_add_ib      $options								>> $log 2>&1
rads_add_orbit   $options -Valt_gdrd --equator --loc-7	>> $log 2>&1
rads_add_ww3_222 $options --all							>> $log 2>&1
rads_add_sla     $options           					>> $log 2>&1

date													>> $log 2>&1

rads_close_sandbox
