
pgrep whereismydata > /dev/null || $(. ./initwhereismydata.sh; date >> $LOG/whereismydata.log ; $PIPE/whereismydata.sh >> $LOG/whereismydata.log 2>&1 ; $PIPE/whereismydata.sh > /ANALYSIS3/monitoring/whereismydata 2>&1)
pgrep delaychart > /dev/null || $(. ./initwhereismydata.sh; $PIPE/delaychart.sh > /ANALYSIS3/log/delaychart 2>&1)
pgrep orbitchart > /dev/null || $(. ./initwhereismydata.sh; $PIPE/orbitchart.sh > /ANALYSIS3/log/orbitchart 2>&1)
#pgrep  -u "$(whoami)" -xf "python ./import_aux.py" > /dev/null || (. ./initpipe36.sh; date >> $LOG/update_monitoring_database.log ; cd $PIPE/import_aux ; python ./import_aux.py >> $LOG/update_monitoring_database.log 2>&1)
#pgrep  -u "$(whoami)" -xf "python ./read_and_import_pl.py" > /dev/null || (. ./initpipe36.sh; date >> $LOG/update_payload_conf.log ; cd $PIPE/payload_conf ; pwd >>  $LOG/update_payload_conf.log ; python ./read_and_import_pl.py >> $LOG/update_payload_conf.log 2>&1)
