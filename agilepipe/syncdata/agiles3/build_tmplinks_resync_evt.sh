#!/bin/bash

datemin='2017-11-15T11:59:59'

tmpevtdir=/tmp/tmpevt_agilepipe
tmpevtdirft=/tmp/tmpevtft_agilepipe

rm -r $tmpevtdir &> /dev/null
mkdir $tmpevtdir
#evt_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='EVT__FM') AND datemin>=\"$datemin\";")

evt_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile f join PIPE_Process p on id=idPipeProcess WHERE (f.type='EVT__FM') AND f.CreatedBy='QLSTDSYSTEM' AND datemin>=\"$datemin\" and not exists (select id from PIPE_Process where status not in ('OK','FAILED') and CreatedBy=p.CreatedBy) ;")

#evt_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='EVT__FM') AND datemin>=DATE_SUB(NOW(), INTERVAL 1 YEAR);")
for file in $evt_files ; do
    ln -s $file $tmpevtdir/$(basename $file)
done

rm -r $tmpevtdirft &> /dev/null
mkdir $tmpevtdirft
#evt_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='EVT__FT3AB') AND datemin>=\"$datemin\";")

evt_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile f join PIPE_Process p on id=idPipeProcess WHERE (f.type='EVT__FT3AB') AND f.CreatedBy='QLSTDSYSTEM' AND datemin>=\"$datemin\" and not exists (select id from PIPE_Process where status not in ('OK','FAILED') and CreatedBy=p.CreatedBy) ;")

#evt_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='EVT__FT3AB') AND datemin>=DATE_SUB(NOW(), INTERVAL 1 YEAR);")
for file in $evt_files ; do
    ln -s $file $tmpevtdirft/$(basename $file)
done

