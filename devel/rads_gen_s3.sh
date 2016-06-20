#!/bin/bash
org=3a.${1}0
fix=3a.${1}1
rm -rf $RADSROOT/data/${org} $RADSROOT/data/${fix}
ls c???/*.nc | rads_gen_s3 -S${org}
cp -pr $RADSROOT/data/${org} $RADSROOT/data/${fix}

# Fix erroneous COG on earlier data
rads_fix_s3 -S${fix} --ymd=20160401000000,20160424230055 --range
# Make the remaining fixes
rads_fix_s3 -S${fix} --all
# Recompute SSB
rads_add_ssb -S${fix} --ssb=ssb_cls
rads_add_ssb -S${fix} --ssb=ssb_cls_c
rads_add_ssb -S${fix} --ssb=ssb_cls_plrm
# Recompute dual freq iono and smooth it
rads_add_dual -S${fix} --recompute
rads_add_dual -S${fix} --recompute --ext=plrm
# Add POE orbit
rads_add_orbit -S${fix} -Valt_cnes --dir=poe
# General geophysical corrections
rads_add_common -S${fix}
rads_add_mog2d -S${fix}
rads_add_ncep -S${fix} -gdwi
rads_add_iono -S${fix} -gn
# Redetermine SSHA
rads_add_sla -S${fix}
rads_add_sla -S${fix} --ext=plrm
