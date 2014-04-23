#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N scrm 
#$ -t 10001-20000
#$ -j y

source parameters

#######################
program=scrm
cmd=${cmd}
#######################

job=${program}_

prefix=${job}${rep}
mkdir ${top_dir}"/"${prefix}

fileprefix=${top_dir}"/"${prefix}"/"${prefix}

{ time -p ${program} ${cmd} > ${fileprefix} -seed ${rep} ;} 2> ${fileprefix}time.text

tree_file_name=${fileprefix}"Trees"
tree_freq_name=${fileprefix}"TreeFreq"
tmrca_name=${fileprefix}"Tmrca"
stat_file_name=${fileprefix}"stat"

grep ';' ${fileprefix} | sed -e "s/\\[.*\\]//g" > ${tree_file_name}
grep ";" ${fileprefix} | sed -e "s/\\[//g" | sed -e "s/\\].*;//g" > ${tree_freq_name}
hybrid-Lambda -gt ${tree_file_name} -tmrca ${tmrca_name}

cat ${fileprefix} | sample_stats > ${stat_file_name}
rm ${fileprefix} ${tree_file_name}

