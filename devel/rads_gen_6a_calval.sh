#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2024  Remko Scharroo
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

# Make a temporary list file
tmplst=`mktemp ${TMPDIR:-/tmp}/rads.XXXXXX`

# Do only latest files when using -d<days>
d0=20000101
case $1 in
	-d*)	days=${1:2}; shift
		d0=`date -u -v -${days}d +%Y%m%d 2>&1` || d0=`date -u --date="${days} days ago" +%Y%m%d` ;;
	-|"")	cat - > "$tmplst"; shift ;;
	*)	[[ ! -d "$1" ]] && cat "$1" > "$tmplst" ;;
esac

# If list not yet made, use find to make it
if [ ! -s "$tmplst" ] && [ $# -gt 0 ] ; then
	# Set bookmark according to $d0
	mrk=$RADSDATAROOT/.bookmark
	TZ=UTC touch -t ${d0}0000 "$mrk"
	# First try only STD files
	find "$@" -name "*STD*.nc" -a -newer "$mrk" | sort > "$tmplst"
	# If still empty, try RED files
	[[ ! -s "$tmplst" ]] && find "$@" -name "*RED*.nc" -a -newer "$mrk" | sort > "$tmplst"
fi

# Get the first file name
dir=`head -n 1 "$tmplst"`

# Exit when no file names are provided
[[ -z $dir ]] && rm -f "$tmplst" && exit

# Determine type (unless already specified)
if test -z $type ; then
	case $dir in
		*/LR*|*P4_2__LR*) type=lr ;;
		*/HR*|*P4_2__HR*) type=hr ;;
	esac
	case $dir in
		*/NR*|*_NR_*) type=${type}nr ;;
		*/ST*|*_ST_*) type=${type}st ;;
		*/NT*|*_NT_*) type=${type}nt ;;
	esac
	case $dir in
		*RMC/*)         type=hrrmc ;;
		*REP/*|*_REP_*) type=${type:0:2}g01 ;;
		*DEV/*|*_DEV_*|*VER/*|*_VER_*) type=${type}d ;;
		*VAL/*|*_VAL_*) type=${type}v ;;
		*)              type=${type}o ;;
	esac
fi

# Process "unadultered" files
rads_open_sandbox "6a.${type}0"

date													>  "$log" 2>&1

# Move list of files
mv -f "$tmplst" "$lst"

# Convert only to RADS, nothing else
rads_gen_s6       $options --min-rec=6 < "$lst"			>> "$log" 2>&1

date													>> "$log" 2>&1

# Now continue with the post-processing
rads_reuse_sandbox "6a.${type}1"

rads_fix_s6       $options --all						>> "$log" 2>&1

# Add MOE orbit (for NRT only)
case $type in
	*nr*) rads_add_orbit  $options -Valt_gdrf --dir=moe_doris	>> "$log" 2>&1 ;;
esac

# For LR, do also _mle3, but only for F09 or earlier
extra=
case $type in
	*lr*) grep -q _G...SEN6 $lst || extra="-x mle3" ;;
esac

# General geophysical corrections
rads_add_common   $options										>> "$log" 2>&1
rads_add_mfwam    $options --all --new							>> "$log" 2>&1
rads_add_iono     $options --all								>> "$log" 2>&1
rads_add_orbit    $options -Valt_gps --dir=jplgpspoe -C5-112	>> "$log" 2>&1
rads_add_orbit    $options -Valt_gps --dir=jplgpsmoe -C113-299	>> "$log" 2>&1
# To support GDR-G with backward compatibility
grep -q _G...SEN6 $lst && rads_add_tide $options --models=fes14		>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options -x -x nr $extra						>> "$log" 2>&1
rads_add_sla      $options -x -x nr $extra						>> "$log" 2>&1

date															>> "$log" 2>&1

rads_close_sandbox
