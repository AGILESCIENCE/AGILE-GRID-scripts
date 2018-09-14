#!/bin/bash
#
# Script per sottomettere un job con root
#
# 

#@ shell = /bin/bash
#@ job_type = serial
#@ environment= COPY_ALL
#@ wall_clock_limit = 240:00:00
#RISORSE PER OGNI TASK
#@ resources = ConsumableCpus(1) 
#ConsumableMemory(5000Mb)
##@ node = 1
##@ tasks_per_node = 1
# in alternativa a task_per_node
##@ total_tasks = 4
#@ error   = job1.$(jobid).err
#@ output  = job1.$(jobid).out
