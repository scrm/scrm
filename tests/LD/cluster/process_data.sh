#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N process_data.sh
#$ -t 1-20
#$ -j y

replicate=10
case="Constantpopsize"

joblist=("ms" \
"scrmwindow100000" "scrmwindow70000" "scrmwindow50000" "scrmwindow30000" \
"scrmwindow10000" "scrmwindow7000" "scrmwindow5000" "scrmwindow3000" \
"scrmwindow1000" "scrmwindow500" "scrmwindow0" \
"macs" "macsretain1000" "macsretain10000" "macsretain30000" \
"macsretain50000" "macsretain70000" "macsretain100000" \
"fastsimcoal")

#rep=$(expr $SGE_TASK_ID - 1)
rep=0
#for rep in $(seq 0 1 16)
    #do
    Job=${case}${joblist[${rep}]}
    JobParamFile=${Job}param
echo -e "case: ${case}\n\
nsam: 6\n\
replicate: ${replicate}\n\
seqlen: 10000001\n\
rho: 4000\n\
job: ${Job}_" > ${JobParamFile}
    python process_data.py ${JobParamFile}
    #done

