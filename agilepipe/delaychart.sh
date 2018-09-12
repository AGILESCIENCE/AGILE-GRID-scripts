#!/bin/bash

# Copyright (c) 2016, AGILE team
# Authors: Andrea Zoli <zoli@iasfbo.inaf.it>,
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

data=/ANALYSIS3/monitoring/delay.dat
outfile=/ANALYSIS3/monitoring/delay.svg
outfile2=/ANALYSIS3/monitoring/delay.png

now=$(date -u +%s)
tmax=$(expr $now)
tmin=$(expr $tmax - 172800)
tmaxiso=$(date -d@$tmax +%FT%T)
tminiso=$(date -d@$tmin +%FT%T)

echo "
set terminal svg font 'Verdana, 8'
set output '$outfile'
set style line 101 lt 1 lc rgb '#BBBBBB' lw 0.5
set grid ls 101
set xlabel 'Time (UTC)' offset 0, -11
set ylabel 'Delay (hours)'
set key out horiz center top width 5

set xdata time
set timefmt '%Y-%m-%dT%H:%M:%S'
set format x '%Y-%m-%d %H:%M:%S'
set xrange ['$tminiso':'$tmaxiso']
set yrange [0:]
set xtics rotate by 90 offset 0,-11
plot '$data' using 1:2 with linespoints pointtype 7 pointsize 0.3 title 'GRID LOG',\
     '$data' using 1:3 with linespoints pointtype 7 pointsize 0.3 title 'GRID EVT',\
     '$data' using 1:4 with linespoints pointtype 7 pointsize 0.3 title 'MCAL 3908'
" | gnuplot

echo "
set terminal pngcairo size 1600,1200 font 'Verdana, 20'
set output '$outfile2'
set style line 1 lc rgb '#8b1a0e' pt 7 ps 0.3 lt 1 lw 3
set style line 2 lc rgb '#5e9c36' pt 7 ps 0.3 lt 1 lw 3
set style line 3 lc rgb '#364b9c' pt 7 ps 0.3 lt 1 lw 3
set style line 101 lt 1 lc rgb '#BBBBBB' lw 1
set grid ls 101
set xlabel 'Time (UTC)' offset 0,-9
set ylabel 'Delay (hours)'
set key out horiz center top width 5

set xdata time
set timefmt '%Y-%m-%dT%H:%M:%S'
set format x '%Y-%m-%d %H:%M:%S'
set xrange ['$tminiso':'$tmaxiso']
set yrange [0:]
set xtics rotate by 90 offset 0,-9
plot '$data' using 1:2 with linespoints title 'GRID LOG' ls 1,\
     '$data' using 1:3 with linespoints title 'GRID EVT' ls 2,\
     '$data' using 1:4 with linespoints title 'MCAL 3908' ls 3
" | gnuplot
