#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N macsretain0
#$ -t 1-1000
#$ -j y

source parameters_preset
program=macs
retain=0
source macs_process.src
