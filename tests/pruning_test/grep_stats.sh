#!/bin/bash

source parameters

#replicate=100
#top_dir="./prune_test"


foo(){
    statfile=${job}_stat
    rm ${statfile}
    tmrca_momentFile=${job}_moment
    rm ${tmrca_momentFile}
    for i in $(seq 1 1 ${replicate})
        do
        prefix=${job}_${i}
        dir=${top_dir}"/"${prefix}
        if [ -d ${dir} ]
            then
            cat ${dir}"/"${prefix}stat >> ${statfile}
            ./compute_moment.py ${dir}"/"${prefix} ${tmrca_momentFile}
        fi
        done

}


job=ms
foo

job=scrm
foo

job=scrmprune50000
foo

job=scrmprune10000
foo
