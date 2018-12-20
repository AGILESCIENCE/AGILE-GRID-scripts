LOG=/ANALYSIS3/log
nohup watch -n 60 sh ./whereismydata.sh > $LOG/watch.txt &
nohup watch -n 60 sh ./piperun.sh > $LOG/watch.txt &
#nohup watch -n 60 sh ./alertmanager.sh > $LOG/watch.txt &
