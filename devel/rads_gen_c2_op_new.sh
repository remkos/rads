#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2021  Remko Scharroo and Eric Leuliette
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

dir=c2.op
rads_open_sandbox $dir

for type in ${types}; do
	mrk=$type/.bookmark
	TZ=UTC touch -t ${d0}0000 "$mrk"
	find $type/ -name "*.nc" -a -newer "$mrk" | sort > "$lst"
	date >>  "$log" 2>&1
	rads_gen_c2_op		$options < "$lst"	>> "$log" 2>&1
done

# General geophysical corrections
rads_add_common   $options					>> "$log" 2>&1
#rads_add_refframe $options					>> "$log" 2>&1
rads_add_iono     $options --all				>> "$log" 2>&1
# Redetermine SSHA
rads_add_sla      $options					>> "$log" 2>&1

date								>> "$log" 2>&1

# Set Navy data aside
files=`egrep 'NOP|IOP' $log | awk '{print $NF}' | awk -F/ '{printf "%s/%s\n",$3,$4}'`
dirs=`echo $files | awk -F/ '{print $1}' | sort | uniq`

pushd $SANDBOX/c2.op/a
for dir in ${dirs[*]} ; do
	mkdir -p $RADSROOT/ext/c2/to_navy/$dir
	for file in ${files[*]} ; do
		ncrename -hOv alt_gdrd,alt_eiggl04s $file $RADSROOT/ext/c2/to_navy/$file >& /dev/null
	done
done
# Remove old Navy data
popd
find to_navy -type f -mtime +30 | xargs rm -f
find to_navy -type d -empty | xargs rmdir

rads_close_sandbox
