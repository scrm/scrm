#!/bin/bash

dir=test-demo
mkdir ${dir}
cd ${dir}
rm *pdf


rep=1000
seqlen=100000

## compare population sturture for a single population data
COMPAREFILE=compareDemo
rm ${COMPAREFILE}

theta=10
r=10


source ../chisq_r.src

source ../ks_r.src

source ../tmrca_r.src

source ../bl_r.src

source ../process_sample_stats.src	


#case 1
echo "case_1" > current_case
rm scrm* scrmInit*
scrm 10 ${rep} -t ${theta} -r ${r} ${seqlen} -eN 0.4 10.01 -eN 1 0.01  -T -L > scrmout
scrmtime
cat scrmout | sample_stats > scrm_stats

scrm 10 ${rep} -t ${theta} -r ${r} ${seqlen} -eN 0.4 10.01 -eN 1 0.01  -T -L -init scrmReadTrees > scrmInitout
scrmInittime
cat scrmInitout | sample_stats > scrmInit_stats

foo

#case 2
echo "case_2" > current_case
rm scrm* scrmInit*
scrm 16 ${rep} -t ${theta} -r ${r} ${seqlen} -G 5.4 -eG 0.4 1 -eN 0.8 15  -T -L > scrmout
scrmtime
cat scrmout | sample_stats > scrm_stats

scrm 16 ${rep} -t ${theta} -r ${r} ${seqlen} -G 5.4 -eG 0.4 1 -eN 0.8 15 -T -L -init scrmReadTrees > scrmInitout
scrmInittime
cat scrmInitout | sample_stats > scrmInit_stats

foo
