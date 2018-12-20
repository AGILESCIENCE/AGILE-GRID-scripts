#pgrep gcndaemon > /dev/null || (. ./initpipe27.sh ; date >> $LOG/gcndaemon_new.log ; nohup $PIPE/gcn/gcndaemon $PIPE/gcn/conf.txt  >> $LOG/gcndaemon_new.log 2>&1 &)

nohup sh update_repository.sh &

nohup sh pipe_manager_submit.sh &

nohup sh pipe_manager_cancel.sh &
