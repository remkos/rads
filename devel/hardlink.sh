#!/bin/bash

# hardlink.sh - Create hardlinks of new files from one directory structure to another

if test $# -le 1 ; then
	cat <<%
hardlink.sh - Create hardlinks of new files from one directory structure to another

syntax: hardlink.sh from_dir to_dir [ -f ] [ find_options ]

where:
	from_dir     : directory where to find new files
	to_dir       : directory where to hardlink new files to (to_dir replaces from_dir in the pathname)
	-f           : make link irrespective of whether the target file is already hardlinked
	find_options : options to add to the find command for example:
	              -mtime -1 : only look for files less than 24 hour old
	              -links 1  : only look for files that are not yet hardlinked

This command first looks for files in from_dir using the find command with find_options.
Each selected file from from_dir is harlinked to the same location in to_dir whenever
the target does not yet exist or is not already hardlinked or when the -f option is used.
Missing directories in to_dir are created if needed.
%
	exit
fi

# Determine if we have gstat, otherwise use stat
stat=$((which gstat 2>/dev/null)||(which stat))

from_dir=$1
to_dir=$2
force=0
shift 2

case "$1" in
	"-f") force=1 ; shift 1 ;;
esac

while read file ; do
	new=${file/#$from_dir/$to_dir}
	if (test $force -eq 1) || (! test -e $new) || (test $($stat --printf '%h\n' $new) -eq 1) ; then
		echo "$file -> $new"
		mkdir -p $(dirname $new)
		ln -f $file $new
	fi
done < <(find $from_dir -type f $*|sort)
