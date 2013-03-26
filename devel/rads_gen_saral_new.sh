#!/bin/bash
#
# Convert SARAL OGDR and IGDR files to RADS
#
# The most recently updated data in the OGDR and IGDR directories
# will be processed, with the IGDRs superceding the OGDRs
#
# Remko Scharroo - 25 Mar 2013
#######################################################################
. radssandbox.sh

rads_open_sandbox sa a
lst=$SANDBOX/rads_gen_saral_new.lst
imrk=igdr/.bookmark
omrk=ogdr/.bookmark

date								>  $log 2>&1

# Process only OGDR data for the last three days (including current)

d0=`date -u -v -7d +%Y%m%d`
TZ=UTC touch -t ${d0}0000 $omrk
find ogdr/c??? -name "SRL_*.nc" -a -newer $omrk | sort > $lst
rads_gen_saral --ymd=$d0 < $lst		>> $log 2>&1

# Now process all IGDR data that came in during the last four days (including current)

d0=`date -u -v -3d +%Y%m%d`
TZ=UTC touch -t ${d0}0000 $imrk
find igdr/c??? -name "SRL_*.nc" -a -newer $imrk | sort > $lst
rads_gen_saral < $lst				>> $log 2>&1

# Do the patches to all data

radsp_iono   $options jpl iri nic	>> $log 2>&1
radsp_common $options				>> $log 2>&1
radsp_ib     $options				>> $log 2>&1
date								>> $log 2>&1

rads_close_sandbox
