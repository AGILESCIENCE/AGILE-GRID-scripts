#!/bin/bash

datemin='2017-11-15T11:59:59'

tmplogdir=/tmp/tmplog_agilepipe
tmpauxdir=/tmp/tmpaux_agilepipe
tmpcordir=/tmp/tmpcor_agilepipe

rm -r $tmplogdir &> /dev/null
mkdir $tmplogdir
log_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='LOG') AND datemin>=\"$datemin\";")
#log_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/qlstd/',filename) FROM PIPE_ProductFile WHERE (type='LOG') AND datemin>=DATE_SUB(NOW(), INTERVAL 1 YEAR);")
for file in $log_files ; do
    ln -s $file $tmplogdir/$(basename $file)
done


#'TCSP', 'OBTUTC'
rm -r $tmpauxdir &> /dev/null
mkdir $tmpauxdir
#aux_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/aux/data/',path,'/',filename) FROM Auxiliary WHERE TYPE in ('EARTH', 'SAS', 'ACSMAN', 'SOE', 'OPM');")
aux_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/aux/data/',path,'/',filename) FROM Auxiliary WHERE TYPE in ('EARTH', 'SAS', 'ACSMAN', 'SOE', 'OPM', 'PMS', 'TCSP', 'TCO') AND status='CHECK_OK' AND ArrivalDate>=DATE_SUB(NOW(), INTERVAL 1 YEAR);")
for file in $aux_files ; do
    ln -s $file $tmpauxdir/$(basename $file)
done
rm -r $tmpcordir &> /dev/null
mkdir $tmpcordir
cor_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/LV1corr','/',filename) FROM IndexL1Corr WHERE datemin>=\"$datemin\";")
for file in $cor_files ; do
    ln -s $file $tmpcordir/$(basename $file)
done

# add 2014 cor files
#cor_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/LV1corr','/',filename) FROM IndexL1Corr WHERE SubType = 08 AND datemin>=\"2014-01-01T00:00:00\" and datemax<=\"2015-01-01T00:00:00\";")
#for file in $cor_files ; do
#    ln -s $file $tmpcordir/$(basename $file)
#done

# add 2009 cor files from 01/05 to 15/05
#cor_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/LV1corr','/',filename) FROM IndexL1Corr WHERE SubType = 08 AND datemin>=\"2009-05-01T00:00:00\" and datemax<=\"2009-05-15T00:0:00\";")
#for file in $cor_files ; do
#    ln -s $file $tmpcordir/$(basename $file)
#done

# add 2012 cor files from 26/11 to 30/11
#cor_files=$(mysql -s -B -hagiles9 -ugs agile2 -e "SELECT CONCAT('/storage1/agile/agile2/LV1corr','/',filename) FROM IndexL1Corr WHERE SubType = 08 AND datemin>=\"2012-11-26T00:00:00\" and datemax<=\"2012-11-30T00:0:00\";")
#for file in $cor_files ; do
#    ln -s $file $tmpcordir/$(basename $file)
#done
