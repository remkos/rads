#-----------------------------------------------------------------------
# Copyright (c) 2011-2024  Remko Scharroo
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
# This bash script provides two functions to aid the processing of GDR data to RADS in a sandbox
#
# rads_open_sandbox - initialise the RADS sandbox
# Usage: rads_open_sandbox sat
#
# rads_reuse_sandbox - move data from RADS sandbox and continue to use it
# Usage: rads_open_sandbox sat [additional-rsync-options]
#
# rads_close_sandbox - move data from RADS sandbox and remove the sandbox
# Usage: rads_close_sandbox [additional-rsync-options]

export RADSDATAROOT=`rads-config --datadir`
export RADSROOT="${RADSDATAROOT%/data*}"

rads_job=`basename $0 .sh`

rads_open_sandbox () {
	rads_sat=$1
	SANDBOX=`mktemp -d ${TMPDIR:-/tmp}/rads.XXXXXX`
	ln -s "$RADSDATAROOT/nml" "$RADSDATAROOT/conf" "$SANDBOX"
	export RADSDATAROOT="$SANDBOX"
	options="-S$rads_sat $opt"
	log=$RADSDATAROOT/${rads_job}-`date -u +%Y%m%d-%H%M%S`.log
	lst=$RADSDATAROOT/${rads_job}-`date -u +%Y%m%d-%H%M%S`.lst
}

rads_reuse_sandbox () {
	rads_new=$1 ; shift
	rsync -aW "$@" "$SANDBOX/$rads_sat/" "$RADSROOT/data/$rads_sat/"
	mkdir -p "$RADSROOT/data/$rads_sat/log"
	gzip < "$log" > "$RADSROOT/data/$rads_sat/log/$(basename $log).gz"
	mv "$SANDBOX/$rads_sat" "$SANDBOX/$rads_new"
	case $rads_sat in
		*.*) ;; # Do nothing when using special directories
		*) (cd "$RADSROOT/tables" ; make $rads_sat web) ;;
	esac
	rads_sat=$rads_new
	options="-S$rads_sat $opt"
}

rads_close_sandbox () {
	egrep -Hi "fault|error" "$log"
	gzip -f "$log"
	export RADSDATAROOT="$RADSROOT/data"
	rsync -aW --remove-source-files "$@" "$SANDBOX/$rads_sat/" "$RADSDATAROOT/$rads_sat/"
	mkdir -p "$RADSDATAROOT/$rads_sat/log"
	mv -f "$log.gz" "$RADSDATAROOT/$rads_sat/log"
	rm -rf "$SANDBOX"
	case $rads_sat in
		*.*) ;; # Do nothing when using special directories
		*) (cd "$RADSROOT/tables" ; make $rads_sat web) ;;
	esac
}
