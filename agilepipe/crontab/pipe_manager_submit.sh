pgrep  -u "$(whoami)" -xf "python ./submit_run.py 1" > /dev/null || (. ./initpipe36.sh ; date >> $LOG/pipe_manager_submit.log ; cd $PIPEMANAGER >> $LOG/pipe_manager_submit.log ; echo $PIPEMANAGER >>  $LOG/pipe_manager_submit.log  ; python ./submit_run.py 1 >> $LOG/pipe_manager_submit.log )

