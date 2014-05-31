#!/bin/bash
dir=test-GROWTH-recomb
mkdir ${dir}
cd ${dir}
rm *pdf


rep=100000

compareGrowth=compareGrowth
rm ${compareGrowth}

theta=100
seqlen=100000
r=10

source ../chisq_r.src

source ../ks_r.src

source ../tmrca_r.src

source ../bl_r.src

source ../process_sample_stats.src	

#case 1 
echo "15_samples_case1" > current_case
rm ms* scrm*

ms 5 ${rep} -t ${theta} -r ${r} ${seqlen} -eG .3 7.0 -eN 0.3 0.5 -T > msout
mstime

scrm 5 ${rep} -t ${theta} -r ${r} ${seqlen} -eG .3 7.0 -eN 0.3 0.5 -T > scrmout
scrmtime


cat msout | sample_stats > ms_stats

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

#cat msout | sample_stats > ms_stats
#cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
#hybrid-Lambda -gt msTrees -tmrca mstmrca
#hybrid-Lambda -gt msTrees -bl msbl


#cat scrmout | sample_stats > scrm_stats
#cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
#hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
#hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 2
echo "5_samples_case2" > current_case
rm ms* scrm*

ms 5 ${rep} -t ${theta} -r ${r} ${seqlen} -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5 -T > msout
mstime

scrm 5 ${rep} -t ${theta} -r ${r} ${seqlen} -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5 -T > scrmout
scrmtime

#cat msout | sample_stats > ms_stats
#cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
#hybrid-Lambda -gt msTrees -tmrca mstmrca
#hybrid-Lambda -gt msTrees -bl msbl


#cat scrmout | sample_stats > scrm_stats
#cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
#hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
#hybrid-Lambda -gt scrmTrees -bl scrmbl
cat msout | sample_stats > ms_stats

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo


##case 3

echo "15_samples_case_3" > current_case
rm ms* scrm*

ms 15 ${rep} -t ${theta} -r ${r} ${seqlen}  -I 2 3 12 -g 1 44.36 -n 2 0.125 -eg 0.03125 1 0.0 -en 0.0625 2 0.05 -ej 0.09375 2 1 -T > msout
mstime

scrm 15 ${rep} -t ${theta} -r ${r} ${seqlen} -I 2 3 12 -g 1 44.36 -n 2 0.125 -eg 0.03125 1 0.0 -en 0.0625 2 0.05 -ej 0.09375 2 1 -T > scrmout
scrmtime

#cat msout | sample_stats > ms_stats
#cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
#hybrid-Lambda -gt msTrees -tmrca mstmrca
#hybrid-Lambda -gt msTrees -bl msbl


#cat scrmout | sample_stats > scrm_stats
#cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
#hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
#hybrid-Lambda -gt scrmTrees -bl scrmbl
cat msout | sample_stats > ms_stats

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

#case 4
echo "15_samples_case_4" > current_case
rm ms* scrm*

ms 15 ${rep} -t ${theta} -r ${r} ${seqlen} -I 2 3 12 -g 1 11.09 -n 1 4.0 -n 2 0.5 -eg 0.125 1 0.0 -en 0.25 2 .2 -ej 0.375 2 1 -T > msout
mstime
scrm 15 ${rep} -t ${theta} -r ${r} ${seqlen} -I 2 3 12 -g 1 11.09 -n 1 4.0 -n 2 0.5 -eg 0.125 1 0.0 -en 0.25 2 .2 -ej 0.375 2 1 -T > scrmout
scrmtime

#cat msout | sample_stats > ms_stats
#cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
#hybrid-Lambda -gt msTrees -tmrca mstmrca
#hybrid-Lambda -gt msTrees -bl msbl


#cat scrmout | sample_stats > scrm_stats
#cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
#hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
#hybrid-Lambda -gt scrmTrees -bl scrmbl

cat msout | sample_stats > ms_stats

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime
foo

