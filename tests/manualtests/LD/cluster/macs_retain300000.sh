#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N macs_retain300000
#$ -t 1-1000
#$ -j y

source parameters_preset
program=macs
retain=300000
source macs_process.src
