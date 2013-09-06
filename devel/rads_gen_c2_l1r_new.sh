#!/bin/bash
#-----------------------------------------------------------------------
# $Id$
#
# Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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
# Convert CryoSat-2 L1R files to RADS
#
# The most recently updated data in the directory SIR_FDM_L1/LATEST
# will be processed.
#
# syntax: rads_gen_c2_l1r_new.sh
#-----------------------------------------------------------------------

. rads_sandbox.sh

rads_open_sandbox c2 a

mrk=$SANDBOX/bookmark
cycle=

# Pick start date: either 4 days ago or the last processed time, whatever is earlier
d0=`date -u -v -4d -v +4H +%Y%m%d`
d2=20`tail -n 1 $RADSROOT/tables/c2a.cyc | cut -c34-39`
[[ $d2 -lt $d0 ]] && d0=$d2

date										>  $log 2>&1

TZ=UTC touch -t ${d0}0000 $mrk
find SIR_FDM_L1/LATEST -name "CS_*.nc" -a -newer $mrk | sort -r | sort -u -t_ -k8,8 > $lst
rads_gen_c2_l1r $options --ymd=$d0 $* < $lst >> $log 2>&1

orbit_opt="-Valt_gdrd --dir=gdr-d-moe"

rads_fix_c2      $options --all				>> $log 2>&1
rads_add_orbit   $options $orbit_opt --equator --loc-7 --rate	>> $log 2>&1
rads_add_orbit   $options -Valt_eig6c		>> $log 2>&1
rads_add_common  $options 					>> $log 2>&1
rads_add_ecmwf   $options --all				>> $log 2>&1
rads_add_ncep    $options -gs               >> $log 2>&1
rads_add_iono    $options --all				>> $log 2>&1
rads_add_mog2d   $options					>> $log 2>&1
rads_add_ww3_222 $options --all				>> $log 2>&1
rads_add_sla     $options                   >> $log 2>&1

date										>> $log 2>&1

# Set Navy data aside
pushd $SANDBOX/c2/a
for dir in c??? ; do
	mkdir -p $RADSROOT/ext/c2/to_navy/$dir
	for file in $dir/*.nc ; do
		ncrename -hOv alt_gdrd,alt_eiggl04s $file $RADSROOT/ext/c2/to_navy/$file >& /dev/null
	done
done
# Remove old Navy data
popd
find to_navy -type f -mtime +30 | xargs rm -f
find to_navy -type d -empty | xargs rmdir

rads_close_sandbox
