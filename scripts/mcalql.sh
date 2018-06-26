#!/usr/bin/env bash

# Copyright (c) 2016, AGILE team
# Authors: Andrea Zoli <zoli@iasfbo.inaf.it>
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

trap exit ERR

if [[ ${#@} != 2 && ${#@} != 3 && ${#@} != 7 ]] ; then
    echo "$0 Bad number of arguments."
    echo "Usage: $0 t0 outputdir [binsize] [dtstart dtstop zoommin zoommax]"
    exit
fi

t0=$1
outdir=$2
binsize=0.1
dtstart=-1000
dtstop=1000
if [[ ${#@} > 2 ]] ; then
    binsize=$3
fi
if [[ ${#@} > 3 ]] ; then
    dtstart=$4
    dtstop=$5
    zoommin=$6
    zoommax=$7
fi

export MCALOUT="$(readlink -m $outdir)/"

mkdir -p ${MCALOUT}{H,RT,misc,grid_offset,count_rate,sci_RM,elcal}

indexfile="/AGILE_PROC3/DATA_ASDC2/INDEX/COR_3908.index"

telemfile=$(python -c "
import sys

indexfile=\"$indexfile\"
t0=$t0

with open(indexfile) as f:
    rows = []
    i = 0
    found = -1
    for line in f:
        cols = line.rstrip('\n').split(' ')
        if t0 > float(cols[1]) and t0 < float(cols[2]):
            found = i
        i += 1
        rows.append(cols)

if found >= 0:
    col = rows[found]
    print col[0]
")

if [[ -z "$telemfile" ]] ; then
	echo "Cannot find $t0 within the 3908 index file. Exiting."
	exit
fi

cd ${MCALOUT}

rm -f trig.list
mcalanalyzer $telemfile -grb trig.list

orbitname=${telemfile##*PKP}
orbitname=${orbitname%%_*}
rootfile="RT${orbitname}_3908.root"

echo "Plotting $rootfile.."

set -x

plotgrb -r ${MCALOUT}/RT/${rootfile} -t0 $t0 $dtstart $dtstop -tbin $binsize -noan -batch -save -z $zoommin $zoommax
