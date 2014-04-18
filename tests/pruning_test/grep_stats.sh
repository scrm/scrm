#!/bin/bash

source parameters

foo(){
    statfile=${job}_stat
    for i in $(seq 1 1 ${replicate})
        do
        prefix=${job}_${i}
        dir=${top_dir}"/"${prefix}
        if [ -d ${dir} ]
            then
            cat ${dir}"/"${prefix}stat >> ${statfile}
	fi
        done

}

job=ms
foo

job=scrm
foo

job=scrmprune
foo
