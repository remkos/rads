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
# Convert latest Jason-3 OGDR and IGDR files to RADS
#
# The most recently updated data in the OGDR and IGDR directories
# will be processed, with the IGDRs superceding the OGDRs
#
# syntax: rads_gen_j3_new.sh
#-----------------------------------------------------------------------
. rads_sandbox.sh

rads_open_sandbox j3

types="ogdr igdr"

date												>  "$log" 2>&1

for type in ${types}; do
    mrk=${type}/.bookmark
    lst=$SANDBOX/rads_gen_j3_tmp_${type}.lst

    case $type in
	    ogdr)
            # Process only OGDR data for the last three days (including current)
            days=2
            d0=$(date -u -v -${days}d +%Y%m%d 2>/dev/null || date -u --date="${days} days ago" +%Y%m%d)
            TZ=UTC touch -t ${d0}0000 $mrk
            find -L ${type}/c[0-8]?? -name "JA3_*.nc" -a -newer $mrk | sort > "$lst"
            ;;
	    igdr)
            # Process all IGDR data that came in during the last four days (including current)
            days=3
            d0=$(date -u -v -${days}d +%Y%m%d 2>/dev/null || date -u --date="${days} days ago" +%Y%m%d)
            TZ=UTC touch -t ${d0}0000 $mrk
            find -L ${type}/c??? -name "JA3_*.nc" -a -newer $mrk | sort > "$lst"
            ;;
    esac
    
    if [ -s "$lst" ]; then
        rads_gen_jason_gdrf --ymd=$d0 $options < "$lst"	>> "$log" 2>&1

    # Add MOE orbit (for OGDR only)
        case $type in
	        ogdr) rads_add_orbit $options -Valt_cnes --dir=gdr-e-moe --equator --rate	>> "$log" 2>&1
		        ;;
        esac

    # Do the patches to all data
        rads_fix_jason    $options --all					>> "$log" 2>&1
        rads_add_common   $options							>> "$log" 2>&1

        if grep -q _2Pg $lst ; then
            # For GDR-G we add the FES2014 model
            rads_add_tide $options --models=fes14           >> "$log" 2>&1
        else
            # For GDR-F we add MLE3 support
            extra="-x mle3 $extra"
        fi

        # Redetermine SSHA
        rads_add_refframe $options -x $extra                >> "$log" 2>&1
        rads_add_sla      $options -x $extra                >> "$log" 2>&1
    fi
    
    date												>> "$log" 2>&1
done

rads_close_sandbox
