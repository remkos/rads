#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2020  Remko Scharroo
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
# Convert Sentinel-6A files to RADS
#
# The data from <directory>/<type>/<cycle(s)> will be processed and put
# into RADS.
# Files will be created in
# $RADSDATAROOT/6a.<type>0 - Unadultered files
# $RADSDATAROOT/6a.<type>1 - RADSified files
#
# syntax: rads_gen_6a_calval.sh <directory>/<types>/<cycles(s)>
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Determine type
dir=$(dirname $1)
case $dir in
*LR/*) type=lr ;;
*HR/*) type=hr ;;
esac
case $dir in
*OPE/*/NR) type=${type}nrt ;;
*OPE/*/ST) type=${type}stc ;;
*OPE/*/NT) type=${type}ntc ;;
*VAL/*/NR) type=${type}nrf ;;
*VAL/*/ST) type=${type}stf ;;
*VAL/*/NT) type=${type}ntf ;;
esac

# Process "unadultered" files
dir=6a.${type}0
rads_open_sandbox "$dir"

date													>  $log 2>&1

find "$@" -name "*RED*.nc" | sort		> "$lst"
rads_gen_s6 	$options --min-rec=6 < "$lst"			>> $log 2>&1
rads_close_sandbox

# Now process do the same again, and do the post-processing
dir=6a.${type}1
rads_open_sandbox $dir

date													>  $log 2>&1

find "$@" -name "*RED*.nc" | sort		> "$lst"
rads_gen_s6 	  $options --min-rec=6 < "$lst"			>> $log 2>&1

# Add MOE orbit (for NRT and STC only)
	case $type in
		*nr*|*st*) rads_add_orbit  $options -Valt_cnes --dir=moe_doris	>> "$log" 2>&1
	esac

# General geophysical corrections
rads_add_common   $options								>> $log 2>&1
rads_add_iono     $options --all						>> $log 2>&1
# Redetermine SSHA
rads_add_sla      $options								>> $log 2>&1

date													>> $log 2>&1

rads_close_sandbox
