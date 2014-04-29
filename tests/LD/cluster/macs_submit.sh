#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N macs 
#$ -t 1-1000
#$ -j y

casefile=toyparam.src
source ${casefile}

#top_dir="toy_test"
top_dir="/well/bsg/joezhu"
rep=$(expr $SGE_TASK_ID )

case=Constantpopsize

program=macs

job=${case}${program}_



#for rep in $(seq 1 1 10)
    #do 
    prefix=${job}${rep}
    mkdir ${top_dir}"/"${prefix}
    fileprefix=${top_dir}"/"${prefix}"/"${prefix}
    cmd="${nsam} ${seqlen} -r 0.0004 -s ${rep}"
    echo ${cmd}    
    { time -p ${program} ${cmd} ;} 2> ${fileprefix}
    
    
    grep "real" ${fileprefix} > ${fileprefix}time.text
    grep "user" ${fileprefix} >> ${fileprefix}time.text
    grep "Tree:" ${fileprefix} | sed "s/,ARG:.*//g" > ${fileprefix}dummy
    
    tree_change_name=${fileprefix}"change"
    tmrca_name=${fileprefix}"Tmrca"
    first_coal_name=${fileprefix}"FirstCoal"

    sed "s/Tree:.*pos://g" ${fileprefix}dummy | sed "s/,len:.*//g" > ${tree_change_name}
    sed "s/Tree:.*TMRCA://g" ${fileprefix}dummy > ${tmrca_name}
    touch ${first_coal_name}    
    ./macs_process.py ${fileprefix} ${seqlen}
    rm ${fileprefix} ${fileprefix}dummy ${tree_change_name}
    #done
