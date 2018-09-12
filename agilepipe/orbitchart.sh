#!/bin/bash

# Copyright (c) 2016, AGILE team
# Authors: Andrea Zoli <zoli@iasfbo.inaf.it>,
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

corpath=/AGILE_PROC3/DATA_ASDC2/COR
tmp=/opt/prod/AGILEPIPE/tmp.dat
if [ "$1" == "" ] ; then
    pmsfile=$(ls /AGILE_PROC3/DATA_ASDC2/AUX/agile_PMS* | tail -n1)
    data_list=/ANALYSIS3/monitoring/orbit_list.dat
    data=/ANALYSIS3/monitoring/orbit.dat
    outfile=/ANALYSIS3/monitoring/orbit.svg
    outfile2=/ANALYSIS3/monitoring/orbit.png

    # write orbit_full_list
    cat $pmsfile | grep -v "#" | grep -v "^$" | sed 's/ *//g' | awk -F'|' '{ system("IFS=\"-\" read -a d <<< \""$2"\" ; year=${d[0]} ; dayofyear=${d[1]} ; hms=${d[2]} ; enddate=$(date -u -d\"$year-01-01 +$dayofyear days -1 day $hms\" +%s) ; ret=$(date -u -d@$enddate \"+%FT%T\") ; echo -n $ret") ; print " "$6" "$8 }' > /ANALYSIS3/monitoring/orbit_full_list.dat
else
    pmsfile=$1
    data=orbit.dat
    data_list=orbit_list.dat
    outfile=orbit.svg
    outfile2=orbit.png
fi

echo "PMS file: $pmsfile"
echo "data list file: $data_list"
echo "data file: $data"
echo "output file: $outfile"

tmin=0
tmax=0
echo -n "" > $data_list
while read line ; do
    endacq=$(echo $line | cut -d'|' -f2)
    IFS='-' read -a endacqv <<< "$endacq"
    year=${endacqv[0]}
    dayofyear=${endacqv[1]}
    enddate=$(date -u -d"$year-01-01 +$dayofyear days -1 day" +%s)
    orbit=$(echo $line | cut -d'|' -f6)
    if [[ $tmin == 0 ]] ; then
        tmin=$enddate
    fi
    tmax=$enddate
    echo "$(date -u -d@$enddate +%F) $orbit" >> $data_list
done <<< "$(echo -e "$(cat "$pmsfile" | grep -v '^#' | grep 'Yes')")"

files=$(ls -v $corpath/*1_0101* | tail -n150 2> /dev/null)
prev=0
counter=0
for file in $files ; do
    obs_id=$(head -n40 "$file" | zcat 2> /dev/null)
    obs_id=${obs_id#*OBS_ID*=*}
    obs_id=${obs_id%%\/*}
    line=$(grep $obs_id $data_list)
    if [ "$line" != "" ] ; then
        sed -i "s/$line/$line acquired/" $data_list
    fi
done

echo -n "" > $data
dates=$(cat $data_list | cut -d' ' -f1 | uniq)
for date in $dates ; do
    acquired=$(cat $data_list | grep acquired | grep $date | wc -l)
    estimate=$(cat $data_list | grep $date | wc -l)
    echo "$date $estimate $acquired" >> $data
done

tmin=$(expr $tmin - 86400)
tmax=$(expr $tmax + 86400)

tmaxiso=$(date -d@$tmax +%F)
tminiso=$(date -d@$tmin +%F)

echo "
set terminal svg font 'Verdana, 8'
set output '$outfile'
set style line 101 lt 1 lc rgb '#BBBBBB' lw 0.5
set xlabel 'Days (UTC date)' offset 0, -7
set ylabel 'Number of orbits (#)'
set nokey
set xdata time
set timefmt '%Y-%m-%d'
set format x '%Y-%m-%d'
set xrange ['$tminiso':'$tmaxiso']
set yrange [0:16]
set xtics rotate by 90 nomirror offset 0,-7
set ytics nomirror
set grid ytics lc rgb '#bbbbbb' lw 1 lt 0

set key horiz center top

set style fill solid 1 noborder
set style histogram rowstacked
set boxwidth 0.7 relative

plot 14 title '' lc rgb '#FF2222',\
     '$data' using 1:2:xtic(1) with boxes lc rgb '#4876b1' title 'Scheduled orbits',\
     '' using 1:3:xtic(1) with boxes lc rgb '#619f3a' title 'Acquired orbits'
     " | gnuplot

echo "
set terminal pngcairo size 1600,1200 font 'Verdana, 20'
set output '$outfile2'
set style line 101 lt 1 lc rgb '#BBBBBB' lw 1
set xlabel 'Days (UTC date)' offset 0,-5
set ylabel 'Number of orbits (#)'
set nokey
set xdata time
set timefmt '%Y-%m-%d'
set format x '%Y-%m-%d'
set xrange ['$tminiso':'$tmaxiso']
set yrange [0:16]
set xtics rotate by 90 nomirror offset 0,-5
set ytics nomirror
set grid ytics ls 101

set key horiz center top

set style fill solid 1 noborder
set style histogram rowstacked
set boxwidth 0.7 relative

plot 14 title '' lc rgb '#FF2222',\
     '$data' using 1:2:xtic(1) with boxes lc rgb '#4876b1' title 'Scheduled orbits',\
     '' using 1:3:xtic(1) with boxes lc rgb '#619f3a' title 'Acquired orbits'
     " | gnuplot
