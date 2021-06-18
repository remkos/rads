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
# This script updates the copyright notice to the year given in the first argument
# E.g.: copyright_update.sh 2016
#
lst=$TMPDIR/copyright_update.lst
grep -l "Copyright (c) 2011-\d\d\d\d  " `find .. -type f | grep -v .git` > $lst
perl -pi -e "s/Copyright \(c\) 2011-\d\d\d\d  /Copyright (c) 2011-$1  /" `cat $lst`
rm -f $lst
