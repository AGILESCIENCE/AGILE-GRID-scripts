#!/bin/bash

#llq -l > joblist_tmp.txt
#cat joblist_tmp.txt | egrep "Job Step Id:|Cmd:" | sed "s/        Job Step Id: //g" | sed "s/               Cmd://" | sed "N;s/\n//"
#rm joblist_tmp.txt

llq -l | egrep "Job Step Id:|Cmd:" | sed "s/        Job Step Id: //g" | sed "s/               Cmd://" | sed "N;s/\n//"
