#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2023  Remko Scharroo
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
# Convert TOPEX/Poseidon Retracked GDR-F files to RADS
#
# syntax: rads_gen_tp_gdrf.sh <directories>
#-----------------------------------------------------------------------
sat=$1
shift

. rads_sandbox.sh

rads_open_sandbox ${sat}.gdrf
lst=$SANDBOX/rads_gen_${sat}_gdrf.lst

date													>  "$log" 2>&1

for tar in "$@"; do
	case "$tar" in
		*.txz) tar -xJf "$tar"; dir=`basename "$tar" .txz` ;;
		*.tgz) tar -xzf "$tar"; dir=`basename "$tar" .tgz` ;;
		*) dir="$tar" ;;
	esac
	ls "$dir"/TP_GPN_*.nc > "$lst"
	rads_gen_tp_gdrf $options < "$lst"					>> "$log" 2>&1
	case "$tar" in
		*.t?z) chmod -R u+w "$dir"; rm -rf "$dir" ;;
	esac
done

# Do the patches to all data
# Moved everything from rads_add_common here, because we need to make some changes

rads_add_grid     $options -Vtopo_srtm15plus					>> "$log" 2>&1
rads_add_era5     $options --all								>> "$log" 2>&1
rads_add_grid     $options -Vdist_coast,gia,basin				>> "$log" 2>&1
rads_add_grid     $options -Vgeoid_egm2008,mss_cnescls15		>> "$log" 2>&1
rads_add_grid     $options -Vmss_dtu15,mss_dtu18				>> "$log" 2>&1
rads_add_grid     $options -Vmss_comb15							>> "$log" 2>&1
rads_add_grid     $options -Vgeoid_eigen6,mss_dtu21				>> "$log" 2>&1
rads_add_grid     $options -Vprox_coast							>> "$log" 2>&1
rads_add_surface  $options										>> "$log" 2>&1
rads_add_surface  $options -s									>> "$log" 2>&1
rads_add_tide     $options --models=stide,ptide,got410,annual	>> "$log" 2>&1
rads_add_tide     $options --models=fes14,lptide				>> "$log" 2>&1
rads_add_tide     $options --models=hret						>> "$log" 2>&1
rads_add_webtide  $options										>> "$log" 2>&1
rads_add_sst      $options --all								>> "$log" 2>&1
rads_add_seaice   $options										>> "$log" 2>&1
rads_add_ncep     $options --dry --wet --air					>> "$log" 2>&1
rads_add_iono     $options --all								>> "$log" 2>&1

case $sat in
tx) extra="-x mle3"
    rads_add_dual $options -C001-235 -r-0.00886 -l				>> "$log" 2>&1
    rads_add_dual $options -C236-480 -r-0.01627 -l				>> "$log" 2>&1
    rads_add_dual $options -C001-235 -r-0.00886 -l $extra		>> "$log" 2>&1
    rads_add_dual $options -C236-480 -r-0.01627 -l $extra		>> "$log" 2>&1
    ;;
esac

rads_add_refframe $options -x $extra							>> "$log" 2>&1
rads_add_sla      $options -x $extra							>> "$log" 2>&1

date															>> "$log" 2>&1

rads_close_sandbox
