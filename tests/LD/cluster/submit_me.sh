#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N redo 
#$ -t 1-10
#$ -j y

case=Constantpopsize
nsam=6
replicate=1000
seqlen=10000001
rho=4000 # This is equal to r * seqlen * 4 * Ne, which is 1e-8 * 1e7 * 4 * 1e4

seqlen=1000001
rho=40 # This is equal to r * seqlen * 4 * Ne, which is 1e-8 * 1e7 * 4 * 1e4
replicate=10

paramfile=${case}param
echo "case: ${case}" >  ${paramfile}
echo "nsam: ${nsam}"  >> ${paramfile}
echo "replicate: ${replicate}" >> ${paramfile}
echo "seqlen: ${seqlen}" >> ${paramfile}
echo "rho: ${rho}" >> ${paramfile}
##echo "diverg

cmd="${nsam} 1 -T -r ${rho} ${seqlen}"

program=ms
job=${case}${program}_
echo "job: ${job}" >> ${paramfile}
rm -r ${job}*

rep=$(expr $SGE_TASK_ID )

#for rep in $(seq 1 1 ${replicate})
    #do
prefix=${job}${rep}
mkdir ${prefix}
fileprefix=${prefix}"/"${prefix}
{ time -p ${program} ${cmd} > ${fileprefix} -seed ${rep} ${rep} ${rep} ;} 2> ${fileprefix}time.text
tree_file_name=${fileprefix}"Trees"
tree_freq_name=${fileprefix}"TreeFreq"
tmrca_name=${fileprefix}"Tmrca"
first_coal_name=${fileprefix}"FirstCoal"

grep ';' ${fileprefix} | sed -e "s/\\[.*\\]//g" > ${tree_file_name}
grep ";" ${fileprefix} | sed -e "s/\\[//g" | sed -e "s/\\].*;//g" > ${tree_freq_name}
hybrid-Lambda -gt ${tree_file_name} -tmrca ${tmrca_name}
hybrid-Lambda -gt ${tree_file_name} -firstcoal ${first_coal_name}
    #done


program=scrm

#exact_window=(1000 0)
exact_window=(100000 50000 10000 1000 0)

for exact_window_i in "${exact_window[@]}"
    do
    job=${case}${program}window${exact_window_i}_
    echo "job: ${job}" >> ${paramfile}
    rm -r ${job}*
    #rm -r ${case}${program}window${exact_window_i}_*
    #for rep in $(seq 1 1 ${replicate})
        #do
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
        #done
    done


./ld_test.py ${paramfile}
