#!/bin/bash

# This script adds the latest ERA5 data to a set of missions.
#
# Syntax: rads-add-era5.sh <days> <mission(s)>
#
# where:
#   <days> is number of days back to consider modifications of the ERA5 data base
#   <mission(s)> are the missions for which to apply the ERA5 corrections (e.g. 6a.{lr,hr}nto1)
#
days=$1
shift 1

# Get the oldest change within the last $days
month=$(cd $ALTIM/data/era5 ; find 2??? -mindepth 1 -maxdepth 1 -type d -mtime -$days | sort | head -n 1 | cut -c1-4,6-7)
[[ ${#month} -eq 0 ]] && exit

cd $RADSROOT/data

for sat in $* ; do
	log=$sat/log/rads_add_era5-$(date +%Y%m%d-%H%M%S).log
	rads_add_era5 -S$sat --ymd=${month}00, --all --new | gzip > $log.gz
done
