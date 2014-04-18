#!/bin/bash

source parameters

foo(){
    statfile=${job}_stat
    for i in $(seq 1 1 ${rep})
        do
        prefix=${job}_${i}
        dir=${top_dir}"/"${prefix}
        if -d dir
            do
            cat ${dir}"/"${prefix}stat >> ${statfile}
            done
        done
}

job=ms
foo

job=scrm
foo

job=scrmprune
foo
