#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N scrm100000 
#$ -t 1-1000
#$ -j y

case=Constantpopsize
nsam=6
replicate=1000
seqlen=10000001
rho=4000 # This is equal to r * seqlen * 4 * Ne, which is 1e-8 * 1e7 * 4 * 1e4

cmd="${nsam} 1 -T -r ${rho} ${seqlen}"
rep=$(expr $SGE_TASK_ID )
#######################
exact_window_i=100000
program=scrm
#######################

job=${case}${program}window${exact_window_i}_
prefix=${case}${program}window${exact_window_i}_${rep}
mkdir ${prefix}
fileprefix=${prefix}"/"${prefix}
{ time -p ${program} ${cmd} > ${fileprefix} -seed ${rep} -l ${exact_window_i} ;} 2> ${fileprefix}time.text
tree_file_name=${fileprefix}"Trees"
tree_freq_name=${fileprefix}"TreeFreq"
tmrca_name=${fileprefix}"Tmrca"
first_coal_name=${fileprefix}"FirstCoal"

grep ';' ${fileprefix} | sed -e "s/\\[.*\\]//g" > ${tree_file_name}
grep ";" ${fileprefix} | sed -e "s/\\[//g" | sed -e "s/\\].*;//g" > ${tree_freq_name}
hybrid-Lambda -gt ${tree_file_name} -tmrca ${tmrca_name}
hybrid-Lambda -gt ${tree_file_name} -firstcoal ${first_coal_name}
