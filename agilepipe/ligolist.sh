#!/bin/bash
#pcregrep -M ".*Got a notice.*\n.*\n.*\n.*LVC*" /ANALYSIS3/log/gcndaemon.log | egrep "notice|LVC" | tail -n80 | paste - - | awk -F' ' '{print $1 $5}' | sed 's|ivorn="ivo://nasa.gsfc.gcn/LVC#| |' | sed 's/"//g'
pcregrep -M ".*Got a notice.*\n.*\n.*\n.*LVC.*\n.*role.*" /ANALYSIS3/log/gcndaemon.log | egrep "notice|LVC|role" | paste - - - | awk -F' ' '{print $1" "$6" "$5}'| sed 's|ivorn="ivo://nasa.gsfc.gcn/LVC#||' | sed 's/"//g'
