#!/bin/bash
cat /ANALYSIS3/log/gcndaemon.log | pcregrep -M ".*Got a notice.*\n.*\n.*\n.*" | egrep "notice|ivorn" | paste - - | awk -F' ' '{print $1 $5}' | sed 's|ivorn="ivo://nasa.gsfc.gcn/| |' |  sed 's/"//g'
