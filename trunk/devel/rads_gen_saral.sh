#!/bin/bash
#-----------------------------------------------------------------------
# $Id: rads_gen_saral.f90 377 2013-03-27 18:58:08Z remko@altimetrics.com $
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
# Convert SARAL O/I/GDR files to RADS
#
# syntax: rads_gen_saral.sh <directories>
#-----------------------------------------------------------------------
. radssandbox.sh

rads_open_sandbox sa a
lst=$SANDBOX/rads_gen_saral_new.lst

date								>  $log 2>&1

for tar in $*; do
	case $tar in
		*.txz) tar -xJf $tar; dir=`basename $tar .txz` ;;
		*.tgz) tar -xzf $tar; dir=`basename $tar .tgz` ;;
		*) dir=$tar ;;
	esac
	ls $dir/SRL_*.nc > $lst
	rads_gen_saral < $lst			>> $log 2>&1
	case $tar in
		*.t?z) chmod -R u+w $dir; rm -rf $dir ;;
	esac
done

# Do the patches to all data

radsp_iono   $options jpl iri nic	>> $log 2>&1
radsp_common $options				>> $log 2>&1
radsp_mog2d  $options				>> $log 2>&1
radsp_ib     $options				>> $log 2>&1

date								>> $log 2>&1

rads_close_sandbox
