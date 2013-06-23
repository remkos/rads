# rads_sandbox.sh
#
# This bash script provides two functions to aid the processing of GDR data to RADS in a sandbox
#
# rads_open_sandbox - initialise the RADS sandbox
# Usage: rads_open_sandbox sat phase
#
# rads_close_sandbox - move data from RADS sandbox and remove the sandbox
# Usage: rads_close_sandbox [additional-rsync-options]

. radsconfig.sh

rads_job=`basename $0 .sh`

rads_open_sandbox () {
rads_sat=$1
rads_phase=$2
SANDBOX=`mktemp -d $RADSROOT/sandbox/XXXXXX`
ln -s $RADSDATAROOT/nml $RADSDATAROOT/conf $SANDBOX
export RADSDATAROOT=$SANDBOX
options="-S$rads_sat/$rads_phase"
log=$RADSDATAROOT/${rads_job}-`date -u +%Y%m%d-%H%M%S`.log
lst=$RADSDATAROOT/${rads_job}-`date -u +%Y%m%d-%H%M%S`.lst
}

rads_close_sandbox () {
egrep -Hi "fault|error" $log
gzip -f $log
export RADSDATAROOT=$RADSROOT/data
rsync -aW --remove-source-files $* $SANDBOX/$rads_sat/ $RADSDATAROOT/$rads_sat/
mv -f $log.gz $RADSDATAROOT/$rads_sat/log
rm -rf $SANDBOX
cd $RADSROOT/tables
make $rads_sat$rads_phase web
}
