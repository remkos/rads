#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2019  Remko Scharroo
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
# This bash script removes from a RADS directory tree all files that
# are older than the file .bookmark at the top of that tree.
#
# This allows to purge all older files after some reprocessing. The
# file .bookmark in that case is 'touched' such that all reprocessed data
# are newer than .bookmark
#
# rads_purge_older.sh -- Remove older files from RADS directory tree
# Usage: rads_purge_older.sh <directory>
#
for dir in $* ; do
	find $dir -not -newer $dir/.bookmark -not -name .bookmark -type f | xargs rm -f
done
