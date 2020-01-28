#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2019  Remko Scharroo
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
# Convert Sentinel-3A SRAL L2 files to RADS
#
# The all standard_measurement.nc files in the named directory will be
# processed.
#
# syntax: rads_gen_3a.sh <directory>
#-----------------------------------------------------------------------
. rads_sandbox.sh

rads_open_sandbox 3a
lst=$SANDBOX/rads_gen_3a.lst

date													>  $log 2>&1

find $* -name "*.nc" | sort > $lst
rads_gen_s3 	  $options --min-rec=6 < $lst			>> $log 2>&1

# General geophysical corrections
rads_add_common   $options								>> $log 2>&1
rads_add_refframe $options --ext=plrm					>> $log 2>&1
rads_add_mog2d    $options								>> $log 2>&1
rads_add_ncep     $options -gdwi						>> $log 2>&1
rads_add_ecmwf    $options -dwui						>> $log 2>&1
rads_add_iono     $options --all						>> $log 2>&1
# Redetermine SSHA
rads_add_sla      $options								>> $log 2>&1
rads_add_sla      $options --ext=plrm					>> $log 2>&1

date													>> $log 2>&1

rads_close_sandbox
