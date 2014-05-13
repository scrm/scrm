#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N ms 
#$ -t 1-1000
#$ -j y

source parameters_preset
#######################
program=ms
job=${case}${program}_
prefix=${job}${rep}
cmd="${cmd} -seed ${rep} ${rep} ${rep}"
#######################
source process_actions.src

