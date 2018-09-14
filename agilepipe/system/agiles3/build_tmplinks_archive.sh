#!/bin/bash

tmplogdir=$HOME/tmplog_agilepipe_archive
tmpevtdir=$HOME/tmpevt_agilepipe_archive

rm -r $tmplogdir &> /dev/null
mkdir $tmplogdir
log_files=$(mysql -s -B -hagiles9 -ugs agile3 -e "SELECT concat('/storage1/agile/agile3/',path,'/',filename) from PIPE_ArchivedFile where cfgsubtype='STD1Kal' and (type='LOG') order by type;")
for file in $log_files ; do
    ln -s /storage1/agile/agile2/qlstd/$file $tmplogdir/$(basename $file)
done
rm -r $tmpevtdir &> /dev/null
mkdir $tmpevtdir
evt_files=$(mysql -s -B -hagiles9 -ugs agile3 -e "SELECT concat('/storage1/agile/agile3/',path,'/',filename) from PIPE_ArchivedFile where cfgsubtype='STD1Kal' and (type='EVT' and subtype='FM') order by type;")
for file in $evt_files ; do
    ln -s /storage1/agile/agile2/qlstd/$file $tmpevtdir/$(basename $file)
done
