#!/bin/bash
#-----------------------------------------------------------------------
# $Id: rads_gen_c2_l1r.sh 759 2015-01-06 13:22:05Z remko@altimetrics.com $
#
# Copyright (c) 2011-2015  Remko Scharroo
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
# This script updates the copyright notice to the year given in the first argument
# E.g.: copyright_update.sh 2015
#
lst=$TMPDIR/copyright_update.lst
for file in `find .. -type f` ; do grep -l "Copyright (c) ....-....  Remko Scharroo (Altimetrics LLC)" $file ; done > $lst
for file in `grep -v copyright_update.sh $lst` ; do
	cp -p $file $file~
	sed "s/Copyright (c) ....-....  Remko Scharroo.*\$/Copyright (c) 2011-$1  Remko Scharroo/" $file~ > $file
done
rm -f $lst
