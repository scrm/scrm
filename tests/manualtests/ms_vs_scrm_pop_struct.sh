#!/bin/bash

dir=test-POP
mkdir ${dir}
cd ${dir}
rm *pdf


rep=100000

## compare population sturture for a single population data
COMPAREFILE=comparePop
rm ${COMPAREFILE}

theta=10


source ../chisq_r.src

source ../ks_r.src

source ../tmrca_r.src

source ../bl_r.src

source ../process_sample_stats.src	


#case 1 
echo "2_samples_case1" > current_case
rm ms* scrm*
ms 2 ${rep} -t ${theta} -eN 0.4 10.01 -eN 1 0.01 -T > msout
scrm 2 ${rep} -t ${theta} -eN 0.4 10.01 -eN 1 0.01  -T > scrmout


cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 1 
echo "2_samples_case1.1" > current_case
rm ms* scrm*
ms 2 ${rep} -t ${theta} -eN 0 1 -eN 0.4 10.01 -eN 1 0.01 -T > msout
scrm 2 ${rep} -t ${theta} -eN 0 1 -eN 0.4 10.01 -eN 1 0.01  -T > scrmout


cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 2
echo "2_samples_case1.2" > current_case
rm ms* scrm*
ms 2 ${rep} -t ${theta} -eN 0 10 -eN 0.4 10.01 -eN 1 0.01 -T > msout
scrm 2 ${rep} -t ${theta} -eN 0 10 -eN 0.4 10.01 -eN 1 0.01  -T > scrmout


cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 2
echo "2_samples_case2" > current_case
rm ms* scrm*
ms 2 ${rep} -t ${theta} -eN .5 0.01 -T > msout
scrm 2 ${rep} -t ${theta} -eN .5 0.01  -T > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo


##case 3

echo "5_samples" > current_case
rm ms* scrm*

ms 5 ${rep} -t ${theta} -eN 0.5 10.0 -T > msout
scrm 5 ${rep} -t ${theta} -eN 0.5 10.0 -T > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 4
echo "4_samples" > current_case
rm ms* scrm*

ms 4 ${rep} -t ${theta} -eN 0.8 15  -T > msout
scrm 4 ${rep} -t ${theta} -eN 0.8 15  -T > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 5
echo "6_samples_case1" > current_case
rm ms* scrm*
ms 6 ${rep} -t ${theta} -eN 1 .1 -eN 3 10 -T > msout
scrm 6 ${rep} -t ${theta} -eN 1 .1 -eN 3 10  -T > scrmout
cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 6
echo "6_samples_case2" > current_case
rm ms* scrm*
ms 6 ${rep} -t ${theta}  -eN .3 10 -T > msout
scrm 6 ${rep} -t ${theta}  -eN .3 10  -T > scrmout
cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo
