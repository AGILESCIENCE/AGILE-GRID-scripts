#!/bin/bash

# Copyright (c) 2016, AGILE team
# Authors: Andrea Zoli <zoli@iasfbo.inaf.it>
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

aux_archive=/AGILE_PROC3/DATA_ASDC2/AUX
log_archive=/AGILE_PROC3/DATA_ASDC2/LOG
evt_archive=/AGILE_PROC3/FM3.119_ASDC2/EVT
cor_archive=/AGILE_PROC3/DATA_ASDC2/COR

export TZ=GMT+0

files=$(ls -v $aux_archive/*.OPM 2> /dev/null)
if [ $? == 0 ] ; then
    file=$(echo "$files" | tail -n1)
    last_time=$(tail -n20 $file | grep Epoch\> | tail -n1 | awk '{print $2, $3}')
    delta_time=$(date -d@$(expr `date +%s -d "$last_time"` - `date +%s`) +"still %d days %H hours")
    timestr=$(time.rb 1 $(echo $last_time | sed 's/\//-/g' | sed 's/ /T/') | head -n3 | tr '\n' ' ')
    printf "OPM   %s -- %-70s -- %s\n" "$delta_time" "$file" "$timestr"
else
    echo OPM   None - -
fi

files=$(ls -v $aux_archive/*.EARTH 2> /dev/null)
if [ $? == 0 ] ; then
    file=$(echo "$files" | tail -n1)
    last_time=$(tail -n20 $file | grep Epoch\> | tail -n1 | awk '{print $2, $3}')
    delta_time=$(date -d@$(expr `date +%s -d "$last_time"` - `date +%s`) +"still %d days %H hours")
    timestr=$(time.rb 1 $(echo $last_time | sed 's/\//-/g' | sed 's/ /T/') | head -n3 | tr '\n' ' ')
    printf "EARTH %s -- %-70s -- %s\n" "$delta_time" "$file" "$timestr"
else
    echo EARTH None - -
fi

files=$(ls -v $aux_archive/*.SAS 2> /dev/null)
if [ $? == 0 ] ; then
    file=$(echo "$files" | tail -n1)
    last_time=$(tail -n20 $file | grep Epoch\> | tail -n1 | awk '{print $2, $3}')
    delta_time=$(date -d@$(expr `date +%s -d "$last_time"` - `date +%s`) +"still %d days %H hours")
    timestr=$(time.rb 1 $(echo $last_time | sed 's/\//-/g' | sed 's/ /T/') | head -n3 | tr '\n' ' ')
    printf "SAS   %s -- %-70s -- %s\n" "$delta_time" "$file" "$timestr"
else
    echo SAS   None - -
fi

files=$(ls -v $aux_archive/*.SOE 2> /dev/null)
if [ $? == 0 ] ; then
    file=$(echo "$files" | tail -n1)
    last_time=$(tail -n20 $file | grep "Nominal_Start\>" | tail -n1 | awk '{print $2, $3}')
    delta_time=$(date -d@$(expr `date +%s -d "$last_time"` - `date +%s`) +"still %d days %H hours")
    timestr=$(time.rb 1 $(echo $last_time | sed 's/\//-/g' | sed 's/ /T/') | head -n3 | tr '\n' ' ')
    printf "SOE   %s -- %-70s -- %s\n" "$delta_time" "$file" "$timestr"
else
    echo SOE   None - -
fi

files=$(ls -v $log_archive/*.LOG.gz 2> /dev/null)
if [ $? == 0 ] ; then
    file=$(echo "$files" | tail -n1)
    last_time=$(head -n20 "$file" | zcat 2> /dev/null)
    last_time=${last_time#*DATE-END*=*\'}
    last_time=${last_time%%\'*}
    last_time=$(echo "$last_time" | sed 's/T/ /' | sed 's/:60/:59/') # :60 for minutes/seconds should not happen, reset if time is bugged..
    delta_time=$(date -d@$(expr `date +%s` - `date +%s -d "$last_time"`) +"%H hours %M min ago   ")
    delta=$(date -d@$(expr `date +%s` - `date +%s -d "$last_time"`) "+%s")
    delta_log=$(python -c "print $delta / 3600.0")
    timestr=$(time.rb 1 $(echo $last_time | sed 's/\//-/g' | sed 's/ /T/') | head -n3 | tr '\n' ' ')
    printf "LOG   %s -- %-70s -- %s\n" "$delta_time" "$file" "$timestr"
else
    echo LOG   None - -
fi

files=$(ls -v $evt_archive/*.EVT__FM.gz 2> /dev/null)
if [ $? == 0 ] ; then
    file=$(echo "$files" | tail -n1)
    last_time=$(head -n20 "$file" | zcat 2> /dev/null)
    last_time=${last_time#*DATE-END*=*\'}
    last_time=${last_time%%\'*}
    last_time=$(echo "$last_time" | sed 's/T/ /' | sed 's/:60/:59/') # :60 for minutes/seconds should not happen, reset if time is bugged..
    delta_time=$(date -d@$(expr `date +%s` - `date +%s -d "$last_time"`) +"%H hours %M min ago   ")
    delta=$(date -d@$(expr `date +%s` - `date +%s -d "$last_time"`) "+%s")
    delta_evt=$(python -c "print $delta / 3600.0")
    timestr=$(time.rb 1 $(echo $last_time | sed 's/\//-/g' | sed 's/ /T/') | head -n3 | tr '\n' ' ')
    printf "EVT   %s -- %-70s -- %s\n" "$delta_time" "$file" "$timestr"
else
    echo EVT   None - -
fi

files=$(ls -v $cor_archive/*_3908_000.lv1.cor.gz 2> /dev/null)
if [ $? == 0 ] ; then
    file=$(echo "$files" | tail -n1)
    last_time=$(head -n20 "$file" | zcat 2> /dev/null)
    last_time=${last_time#*DATE-END*=*\'}
    last_time=${last_time%%\'*}
    last_time=$(echo "$last_time" | sed 's/T/ /' | sed 's/:60/:59/') # :60 for minutes/seconds should not happen, reset if time is bugged..
    delta_time=$(date -d@$(expr `date +%s` - `date +%s -d "$last_time"`) +"%H hours %M min ago   ")
    delta=$(date -d@$(expr `date +%s` - `date +%s -d "$last_time"`) "+%s")
    delta_3908=$(python -c "print $delta / 3600.0")
    timestr=$(time.rb 1 $(echo $last_time | sed 's/\//-/g' | sed 's/ /T/') | head -n3 | tr '\n' ' ')
    printf "3908   %s -- %-70s -- %s\n" "$delta_time" "$file" "$timestr"
else
    echo 3908   None - -
fi

mkdir -p $PATH_RES/monitoring
echo "$(date -u +%FT%T) $delta_log $delta_evt $delta_3908" >> $PATH_RES/monitoring/delay.dat
