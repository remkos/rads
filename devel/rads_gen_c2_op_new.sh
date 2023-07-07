#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2022  Remko Scharroo and Eric Leuliette
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
# Convert latest CryoSat-2 NOP, IOP, and GOP files to RADS
#
# The most recently updated data in the nop, iop, and gop directories
# will be processed.
#
# syntax: rads_gen_c2_op_new.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

# Process NOP/IOP/GOP data
days=${days:-3}
types="${types:-nop iop gop}"
while getopts "nigd:" arg; do
	case $arg in
		d) days=$OPTARG ;;
		n) types=nop ;;
		i) types=iop ;;
		g) types=gop ;;
	esac
done

d0=`date -u -v -${days}d +%Y%m%d 2>&1` || d0=`date -u --date="${days} days ago" +%Y%m%d`

dir=c2
rads_open_sandbox $dir

for type in ${types}; do
	mrk=$type/.bookmark
	TZ=UTC touch -t ${d0}0000 "$mrk"
	find $type/ -name "*.nc" -a -newer "$mrk" | sort > "$lst"
	date >>  "$log" 2>&1
	if [ "$type" == "gop" ]; then
		rads_gen_c2_op 	$options < "$lst"	>> "$log" 2>&1
	else
		rads_gen_c2_op --ymd=$d0 	$options < "$lst"	>> "$log" 2>&1
	fi
done

# General geophysical corrections
rads_add_common   $options					>> "$log" 2>&1
# Include GOT4.8 model for Navy
rads_add_tide     $options --models=got48	>> "$log" 2>&1
rads_add_iono     $options --all			>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options					>> "$log" 2>&1
rads_add_sla      $options					>> "$log" 2>&1

date										>> "$log" 2>&1

# Set Navy and NHC data aside
files_iop=`grep IOP $log | grep written | awk '{print $NF}' | awk -F/ '{printf "%s/%s\n",$3,$4}' | sort | uniq`
files_nop=`grep NOP $log | grep written | awk '{print $NF}' | awk -F/ '{printf "%s/%s\n",$3,$4}' | sort | uniq`

pushd $SANDBOX/c2/a
for file in ${files_iop[*]} ; do
	bfile=`basename $file`
	cp $file $RADSROOT/ext/c2/to_navy/igdr/$bfile
done
for file in ${files_nop[*]} ; do
	bfile=`basename $file`
	cp $file $RADSROOT/ext/c2/to_navy/ogdr/$bfile
done

# Remove old Navy/NHC data
popd
find to_navy -type f -mtime +7 | xargs rm -f

rads_close_sandbox
