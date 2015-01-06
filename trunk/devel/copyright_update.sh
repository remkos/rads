#!/bin/sh
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
