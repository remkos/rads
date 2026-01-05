#!/bin/bash
#-----------------------------------------------------------------------
# Copyright (c) 2011-2026  Remko Scharroo and Eric Leuliette
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
# ogdrsplit.sh -- Split Jason-3 OGDR-F files into pass files
#
# Splits and merges Jason-3 OGDR-F files into separate pass files of
# 1-Hz data. Works for either GDR-D or GDR-F standards.
#
# The pass files will be named <destdir>/cCCC/???_*_2P?PCCC_[P]PPP.nc, where
# CCC is the cycle number and [P]PPP the pass number. The directory
# <destdir>/cCCC will be created if needed.

# The input file names are read from standard input.
# syntax: ogdrsplit.sh <destdir> < <OGDR filenames>
#
# This program needs an up-to-date files $RADSROOT/ext/??/JA?_ORF.txt with
# equator crossing information.
#-----------------------------------------------------------------------

destdir=$1
ogdrs=( )
while (( $# )); do
    ogdrs+=( "$2" ); shift
done

#declare -a lon lat cycle pass eqlon eqtime index starttime startindex

sec1985to2000=473299200
sec1970to2000=946684800
tmp1=$(mktemp -t ogdrftime1)
tmp2=$(mktemp -t ogdrftime2)
tmpdir=$(mktemp -d -t ogdrf)

#Hard coded for Jason-3 for now
filetype='JA3_OPN_2Pf'
sat='JA3'
nr_passes=254
abs_pass_offset=0
case $sat in
JA1)
    orf=${RADSROOT}/ext/j1/JA1_ORF.txt
;;
JA2)
    orf=${RADSROOT}/ext/j2/JA2_ORF.txt
;;
JA3)
    orf=${RADSROOT}/ext/j3/JA3_ORF.txt
;;
*)
  echo "Error"
;;
esac

#Read ORF file
hash=0
iline=0
while read -r ymd hms ocycle opass orbitnr ilon ilat
do
    if [[ $hash -lt 5 ]]; then
        if [[ ${ymd::1} == "#" ]] ; then ((hash+=1)); fi
        continue
    fi
    echo $ymd $hms >> $tmp1
    ocycle=$((10#$ocycle))
    opass=$((10#$opass))
    cycle[$iline]=$ocycle
    pass[$iline]=$opass
    secdec[$iline]=${hms#*.}
    lon[$iline]=$ilon
    lat[$iline]=$ilat
    timeposix=$ymd' '$hms
    timeposix=${timeposix//\//-}
    timeposix[$iline]=$timeposix
    abs_pass[$iline]=$(($(($(($((ocycle - 1)) * nr_passes)) + opass)) + abs_pass_offset))
    ((iline+=1))
done < $orf

#Convert ORF times to epoch 1970 seconds
gdate -u +%s --file=$tmp1 > $tmp2

#Read times and find equator crossings and turning points
iline=0
npass=0
while read -r itime
do
#  sec2000=$((itime - sec1970to2000))\.${secdec[$iline]}
  ilat=${lat[$iline]}
  if [[ "$ilat" =~ ^-?[0](\.0+)?$ ]]; then
      index[$npass]=$iline
      ((npass+=1))
  else
      startindex[$npass]=$iline
      starttime[$npass]=$itime
  fi
  ((iline+=1))
done < $tmp2
rm $tmp2

for ogdr in ${ogdrs[@]}; do
#Find times of first and last data in OGDR
#mapfile only in bash4
#mapfile -t ogdr_times < <( ncks -A -x --cdl $ogdr | egrep 'first|last' )
    ogdr_times=()
    while IFS= read -r line; do
        ogdr_times+=("$line")
    done < <(ncks -A -x --cdl $ogdr | egrep 'first|last')
    ogdr_first_meas_time=${ogdr_times[0]}
    ogdr_first_meas_time=${ogdr_first_meas_time#*'"'}; ogdr_first_meas_time=${ogdr_first_meas_time%'"'*}
    ogdr_first_meas_time_sec="$(gdate +%s --date="$ogdr_first_meas_time UTC")"

    ogdr_last_meas_time=${ogdr_times[1]}
    ogdr_last_meas_time=${ogdr_last_meas_time#*'"'}; ogdr_last_meas_time=${ogdr_last_meas_time%'"'*}
    ogdr_last_meas_time_sec="$(gdate +%s --date="$ogdr_last_meas_time UTC")"

#Search identify passes in the OGDR from turning point times
#Should be equator crosses?
    for i in ${!starttime[@]}; do
        [[ $ogdr_first_meas_time_sec -ge ${starttime[$i]} && $ogdr_first_meas_time_sec -le ${starttime[(($i+1))]} ]] && firstpass=$i
        [[ $ogdr_last_meas_time_sec -ge ${starttime[$i]} && $ogdr_last_meas_time_sec -le ${starttime[(($i+1))]} ]] &&  lastpass=$i
    done

    for ((ipass=${firstpass};ipass<=${lastpass};ipass++)); do
        ieq=${index[$ipass]}
        ist=${startindex[$ipass]}
        ist2=${startindex[(($ipass+2))]}

        printf -v outdir '%s/c%03d' "${destdir}" "${cycle[$ieq]}"
        mkdir -p $outdir
        printf -v outfn '%sP%03d_%03d.nc' "$filetype" "${cycle[$ieq]}" "${pass[$ieq]}"
        outnm=$outdir/$outfn

#By default split with ORF turning point times
        orfdate1=${starttime[$ipass]}
        orfdate2=${starttime[(($ipass+1))]}
#Append decimals for ncks to interpret as a value instead of an index.
        time1=$((orfdate1 - sec1970to2000))'.'${secdec[$ist]}
        time2=$((orfdate2 - sec1970to2000))'.'${secdec[$ist2]}

        exists=0
        if [[ -e "$outnm" ]]; then
            echo "$outnm exists."
            exists=1
            file_last_meas_time=$(ncks -A -x --cdl $outnm | egrep 'last')
            file_last_meas_time=${file_last_meas_time#*'"'}; file_last_meas_time=${file_last_meas_time%'"'*}
            file_secdec=${file_last_meas_time#*.}
            file_last_meas_time_sec="$(gdate +%s --date="$file_last_meas_time UTC")"

            if [[ $ogdr_last_meas_time_sec > $file_last_meas_time_sec && $orfdate2 > $file_last_meas_time_sec ]]; then
                time1=$(($((file_last_meas_time_sec - sec1970to2000)) + 1))'.'$file_secdec
            else
                continue
            fi
        fi

        echo $outnm
        tmpnm=$tmpdir/$outfn

# Remove data_20 group and slice by time range
# Could test for data_20 to make generic for SSHA, SARAL, or GDR-D
        ncks -O --cnk_csh=1000000000 -x -g data_20 -h -d time,$time1,$time2 $ogdr $tmpnm
# Include a record dimension for NCO tools to catenate along
        ncks -h -O --mk_rec_dmn time $tmpnm $tmpnm
        ncatted -O -h -a cycle_number,global,o,i,${cycle[$ieq]} $tmpnm
        ncatted -O -h -a pass_number,global,o,i,${pass[$ieq]} $tmpnm
        ncatted -O -h -a absolute_pass_number,global,o,i,${abs_pass[$ieq]} $tmpnm
        ncatted -O -h -a equator_lon,global,o,d,${lon[$ieq]} $tmpnm
        ncatted -O -h -a equator_time,global,o,c,"${timeposix[$ieq]}" $tmpnm
        ncatted -O -h -a first_meas_time,global,o,c,"${timeposix[$ist]}" $tmpnm
        if [[ "$ipass" == $lastpass ]]; then
            ncatted -O -h -a last_meas_time,global,o,c,"$ogdr_last_meas_time" $tmpnm
        else
            ncatted -O -h -a last_meas_time,global,o,c,"${timeposix[(($ist+2))]}" $tmpnm
        fi
        ncatted -O -h -a original,global,o,c,$(basename $ogdr) $tmpnm

        if [[ "$exists" == "1" ]]; then
# Appending with ncks rather than ncrcat is more efficient
            ncks -A --cnk_csh=1000000000 -h $tmpnm $outnm
#	ncrcat --cnk_csh=1000000000 -O -h $outnm $tmpnm $tmpnm 2>/dev/null
        fi
        \mv -f $tmpnm $outnm
    done
done
#Clean up
\rm -rf $tmp1 $tmpdir
