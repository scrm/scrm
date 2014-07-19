#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N macs_retain30000
#$ -t 1-1000
#$ -j y

source parameters_preset
program=macs
retain=30000
source macs_process.src
