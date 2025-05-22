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
# Convert Sentinel-3 SRAL L2 files to RADS
#
# The all standard_measurement.nc files in the named directory will be
# processed.
#
# syntax: rads_gen_3.sh <directory>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Which satellite?
case $0 in
	*rads_gen_3a*) sat=3a ;;
	*rads_gen_3b*) sat=3b ;;
	*) echo "$0: unknown script" ; exit ;;
esac

# Exit when no directory names are provided
[[ $# -eq 0 ]] && exit

rads_open_sandbox ${sat}

date													>  "$log" 2>&1

find "$@" -name "*.nc" | sort > $lst
rads_gen_s3 	  $options --min-rec=6 < $lst			>> "$log" 2>&1

# Add MOE orbit (for NRT and STC only) and CPOD POE (for NTC/REP only)
case $type in
	nr*|st*) rads_add_orbit  $options -Valt_cnes --dir=moe_doris	>> "$log" 2>&1 ;;
	*)	     rads_add_orbit  $options -Valt_cpod					>> "$log" 2>&1 ;;
esac

# General geophysical corrections
rads_add_common   $options								>> "$log" 2>&1
case ${sat} in
	3a) rads_add_mfwam $options -C40-199 --all --new	>> "$log" 2>&1 ;;
	3b) rads_add_mfwam $options -C21-199 --all --new	>> "$log" 2>&1 ;;
esac
# Redetermine SSHA
rads_add_refframe $options -x -x plrm					>> "$log" 2>&1
rads_add_sla      $options -x -x plrm					>> "$log" 2>&1

date													>> "$log" 2>&1

rads_close_sandbox
