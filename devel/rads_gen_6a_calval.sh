#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2021  Remko Scharroo
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

# Do only latest files when using -d<days>
d0=20000101
case $1 in
	-d*) days=${1:2}; shift
		d0=`date -u -v -${days}d +%Y%m%d 2>&1` || d0=`date -u --date="${days} days ago" +%Y%m%d`
		;;
esac

# Exit when no directory names are provided
[[ $# -eq 0 ]] && exit

# Determine type
dir=$(dirname $1)
red=STD
case $dir in
	*LR/*|*P4_2__LR*) type=lr ;;
	*HR/*|*P4_2__HR*) type=hr ;;
esac
case $dir in
	*/NR*) type=${type}nr; red=RED ;;
	*/ST*) type=${type}st ;;
	*/NT*) type=${type}nt ;;
esac
case $dir in
	*OPE/*) type=${type}o ;;
	*VAL/*) type=${type}v ;;
	*DEV/*) type=${type}d ;;
	*REP/*) type=${type:0:2}rep ;;
	*RMC/*) type=hrrmc ;;
esac

# Process "unadultered" files
rads_open_sandbox "6a.${type}0"

date													>  "$log" 2>&1

# Set bookmark according to $d0
mrk=$RADSDATAROOT/.bookmark
TZ=UTC touch -t ${d0}0000 "$mrk"
find "$@" -name "*${red}*.nc" -a -newer "$mrk" | sort	>  "$lst"

# Convert only to RADS, nothing else
rads_gen_s6       $options --cal1 --min-rec=6 < "$lst"	>> "$log" 2>&1

date													>> "$log" 2>&1

# Now continue with the post-processing
rads_reuse_sandbox "6a.${type}1"

rads_fix_s6       $options --all						>> "$log" 2>&1

# Add MOE orbit (for NRT only)
case $type in
	*nr*) rads_add_orbit  $options -Valt_gdrf --dir=moe_doris	>> "$log" 2>&1 ;;
esac

# For LR, add mle3
case $type in
	*lr*) extra="-x mle3" ;;
	   *) extra= ;;
esac

# General geophysical corrections
rads_add_common   $options								>> "$log" 2>&1
rads_add_mfwam    $options --all						>> "$log" 2>&1
rads_add_iono     $options --all						>> "$log" 2>&1
rads_add_mog2d    $options								>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options -x $extra					>> "$log" 2>&1
rads_add_sla      $options -x $extra					>> "$log" 2>&1

date													>> "$log" 2>&1

rads_close_sandbox
