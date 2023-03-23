#!/bin/bash

# hardlink.sh - Create hardlinks of new files from one directory structure to another

if test $# -le 1 ; then
	cat <<%
$0 - Create hardlinks of new files from one directory structure to another

syntax: $0 from_dir to_dir [ find_options ]

where:
	from_dir     : directory where to find new files
	to_dir       : directory where to hardlink new files to (to_dir replaces from_dir in the pathname)
	find_options : options to add to the find command for example:
	               -mtime -1 : only look for files less than 24 hour old
	               -links 1  : only look for files that are not yet hardlinked

This command first looks for files in from_dir using the find commond with find_options.
Each file that does not exist in the same place in to_dir, or has the same or old time step in to_dir
is then hardlinked from from_dir to to_dir.
Missing directories are created if needed.
%
	exit
fi

from_dir=$1
to_dir=$2
shift 2

if test -x /opt/local/bin/gstat ; then
	stat=/opt/local/bin/gstat
else
	stat=stat
fi

while read file ; do
	new=${file/#$from_dir/$to_dir}
	if ! test $file -ot $new ; then
		echo "$file -> $new"
		mkdir -p $(dirname $new)
		ln -f $file $new
	fi
done < <(find $from_dir -type f $*)
