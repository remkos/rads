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
# Convert REAPER files to RADS
#
# syntax: rads_gen_reaper.sh <directories>
#-----------------------------------------------------------------------
case $1 in
	*ers1*) ers=1 ;;
	*ers2*) ers=2 ;;
	*) echo "Unknown satellite" ; exit ;;
esac
case $1 in
	*RP01*) prod=RP01 ;;
	*) echo "Unknown product" ; exit ;;
esac

. rads_sandbox.sh

rads_open_sandbox e$ers
lst=$SANDBOX/rads_gen_reaper.lst

date											>  "$log" 2>&1

find $1 -name "E${ers}_REAP_ERS_ALT_2M_*_${prod}.nc" > "$lst"

rads_gen_reaper  < "$lst"							>  "$log" 2>&1

# Do the patches to all data

rads_fix_reaper   $options --tide --tbias			>> "$log" 2>&1
rads_add_flags    $options --file=$RADSROOT/ext/reaper/e${1}_flags.dat	>> "$log" 2>&1
rads_add_orbit    $options -Valt_reaper_deos		>> "$log" 2>&1
rads_add_ssb      $options --ssb=ssb_hyb			>> "$log" 2>&1
case $ers in
1) rads_add_iono  $options --nic09					>> "$log" 2>&1 ;;
2) rads_add_iono  $options -C0-34 --nic09			>> "$log" 2>&1
   rads_add_iono  $options -C35-85 --nic09 --gim 	>> "$log" 2>&1 ;;
esac
rads_add_mog2d    $options							>> "$log" 2>&1
rads_add_common   $options							>> "$log" 2>&1
# Redetermine SSHA
rads_add_refframe $options							>> "$log" 2>&1
rads_add_sla      $options							>> "$log" 2>&1

date												>> "$log" 2>&1

rads_close_sandbox
