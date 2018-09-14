#!/bin/bash

datemin='2016-05-15T11:59:59'

tmplogdir=$HOME/tmplog_agilepipe
tmpevtdir=$HOME/tmpevt_agilepipe
tmpauxdir=$HOME/tmpaux_agilepipe
tmpcordir=$HOME/tmpcor_agilepipe

rm -r $tmplogdir &> /dev/null
mkdir $tmplogdir
log_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='LOG') AND datemin>=\"$datemin\";")
#log_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='LOG') AND datemin>=DATE_SUB(NOW(), INTERVAL 1 YEAR);")
for file in $log_files ; do
    ln -s $file $tmplogdir/$(basename $file)
done
rm -r $tmpevtdir &> /dev/null
mkdir $tmpevtdir
evt_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='EVT__FM') AND datemin>=\"$datemin\";")
#evt_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='EVT__FM') AND datemin>=DATE_SUB(NOW(), INTERVAL 1 YEAR);")
for file in $evt_files ; do
    ln -s $file $tmpevtdir/$(basename $file)
done
rm -r $tmpauxdir &> /dev/null
mkdir $tmpauxdir
#aux_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/aux/data/',path,'/',filename) FROM Auxiliary WHERE TYPE in ('EARTH', 'SAS', 'ACSMAN', 'SOE', 'OPM');")
aux_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/aux/data/',path,'/',filename) FROM Auxiliary WHERE TYPE in ('EARTH', 'SAS', 'ACSMAN', 'SOE', 'OPM', 'PMS') AND ArrivalDate>=DATE_SUB(NOW(), INTERVAL 1 YEAR);")
for file in $aux_files ; do
    ln -s $file $tmpauxdir/$(basename $file)
done
rm -r $tmpcordir &> /dev/null
mkdir $tmpcordir
cor_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/LV1corr','/',filename) FROM IndexL1Corr WHERE datemin>=\"$datemin\";")
for file in $cor_files ; do
    ln -s $file $tmpcordir/$(basename $file)
done
