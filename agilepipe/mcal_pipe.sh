#!/usr/bin/env bash

# Copyright (c) 2016, AGILE team
# Authors: Andrea Zoli <zoli@iasfbo.inaf.it>
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

agilepipe_path="$AGILE/AGILEPIPE/"

trap exit ERR

if [[ ${#@} != 1 ]] ; then
    echo "$0 Bad number of arguments."
    echo "Usage: $0 file.conf"
    exit
fi

unset conf
while IFS= read -r line; do
    conf+=("$line")
done < $1

analysistype=${conf[0]}
instrumentname=${conf[1]}
trigid=${conf[2]}
seqnum=${conf[3]}
queue=${conf[4]}
outdir=${conf[5]}
outname=${conf[6]}
t0=${conf[7]}
binsize=${conf[8]}
dtstart=${conf[9]}
dtstop=${conf[10]}
zoommin=${conf[11]}
zoommax=${conf[12]}

rundir="${PATH_RES}/${instrumentname}_RUN/$outname"
mkdir -p $rundir
outdir="${PATH_RES}/${instrumentname}/$outdir"
mkdir -p $outdir

$agilepipe_path/mcalql.sh $t0 $rundir $binsize $dtstart $dtstop $zoommin $zoommax

img=$(ls -1 $rundir/*.gif | head -n1)
triglist=$rundir/trig.list
cp $rundir/img1.gif $outdir/${outname}.gif
cp $triglist $outdir/${outname}_trig.list
