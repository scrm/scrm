#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N scrm700000 
#$ -t 1-1000
#$ -j y


source parameters_preset
#######################
exact_window_i=700000
program=scrm
job=${case}${program}window${exact_window_i}_
prefix=${case}${program}window${exact_window_i}_${rep}
cmd="${cmd} -l ${exact_window_i} -seed ${rep} "
#######################
source process_actions.src

