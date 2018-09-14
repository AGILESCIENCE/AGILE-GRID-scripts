#!/bin/sh


LOG=/ANALYSIS3/log

. ~/.bashrc
cd $AGILE/DeepVar/import/dwh
ruby ./import_fact_table.rb 128 >> $LOG/import_dwh_128.log 2>&1
