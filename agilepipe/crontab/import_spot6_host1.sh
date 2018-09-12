#!/bin/sh



LOG=/ANALYSIS3/log

. ~/.bashrc
cd $AGILE/DeepVar/import

ruby ./import_last_spot6_host1.rb >> $LOG/import_detection_spot6_host1.log 2>&1
