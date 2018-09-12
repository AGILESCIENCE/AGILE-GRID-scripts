#!/bin/bash

# Copyright (c) 2016, AGILE team
# Authors: Andrea Zoli <zoli@iasfbo.inaf.it>,
#          Andrea Bulgarelli <zoli@iasfbo.inaf.it>
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

trap exit SIGINT SIGTERM

TIMEOUT=5
MAX_RETRIES=12

retry()
{
    i=0
    false
    while [[ $? -ne 0 && $i -lt $MAX_RETRIES ]] ;
    do
        i=$(($i+1))
        $@
    done
    if [ $i -eq $MAX_RETRIES ]
    then
        date
        echo "Err: reached the maximum number of $MAX_RETRIES retries with $TIMEOUT seconds of timeout!"
        exit 1
    fi
}

ssh_cmd()
{
    ssh -o ConnectTimeout=$TIMEOUT -o ConnectionAttempts=$MAX_RETRIES $@
    if [[ $? -ne 0 ]] ; then
        date
        echo "Err: reached the maximum number of $MAX_RETRIES retries with $TIMEOUT seconds of timeout!"
        exit 1
    fi
}

echo "Rebuild gs@agiles3 symbolic links.."
ssh_cmd gs@agiles3 '$HOME/build_tmplinks_archive.sh'
echo "Rebuild complete."

date
echo "Copying files from gs@agiles3.."
retry rsync -rLptgoDvz --delete --timeout=$TIMEOUT 'gs@agiles3:$HOME/tmplog_agilepipe_archive/' '/ASDC_PROC3/DATA_ASDCSTDk/LOG/'
retry rsync -rLptgoDvz --delete --timeout=$TIMEOUT 'gs@agiles3:$HOME/tmpevt_agilepipe_archive/' '/ASDC_PROC3/FM3.119_ASDCSTDk/EVT/'
echo "Transfer complete."

echo "Rebuild the indexes.."
log_index=/ASDC_PROC3/DATA_ASDCSTDk/INDEX/LOG.log.index
AG_indexgen /ASDC_PROC3/DATA_ASDCSTDk/LOG LOG ${log_index}_tmp
sort -k3 ${log_index}_tmp -o ${log_index}_tmp
mv ${log_index}_tmp ${log_index}
evt_index=/ASDC_PROC3/FM3.119_ASDCSTDk/INDEX/EVT.index
AG_indexgen /ASDC_PROC3/FM3.119_ASDCSTDk/EVT EVT ${evt_index}_tmp
sort -k3 ${evt_index}_tmp -o ${evt_index}_tmp
mv ${evt_index}_tmp ${evt_index}
echo "Rebuild complete."

date
