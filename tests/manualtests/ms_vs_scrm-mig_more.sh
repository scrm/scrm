#!/bin/bash

#Compare summary statistics of ms and scrm for migration

dir=test-mig-more
mkdir ${dir}
cd ${dir}
rm *pdf


rep=100000

## compare TMRCA
COMPAREFILE=compareMIG
rm ${COMPAREFILE}


theta=10

source ../chisq_r.src

source ../ks_r.src

source ../tmrca_r.src

source ../bl_r.src

source ../process_sample_stats.src	
	
	
#case 1 
#2 sub population, 1 sample from each subpopulation, mutation rate is 5
echo "2groups1sam1sam_case1_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 2 ${rep} -t ${theta} -I 2 1 1 5.0 -T -L > msout
scrm 2 ${rep} -t ${theta} -I 2 1 1 -ma x 5.0 5.0 x -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

#case 2
#2 sub population, 1 sample from each subpopulation, mutation rate is 5
echo "2groups1sam1sam_case2_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 2 ${rep} -t ${theta} -I 2 1 1 5.0 -T -L > msout
scrm 2 ${rep} -t ${theta} -I 2 1 1 5.0 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo


#case 3
#2 sub population, 1 sample from each subpopulation, mutation rate is 5
echo "2groups1sam1sam_case3_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 2 ${rep} -t ${theta} -I 2 1 1 -ma x 5 0 x -ej 2.0 2 1 -T -L > msout
scrm 2 ${rep} -t ${theta} -I 2 1 1 -ma x 5 0 x -ej 2.0 2 1 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

#case 3
#2 sub population, 1 sample from each subpopulation, mutation rate is 5
echo "2groups1sam1sam_case4_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 2 ${rep} -t ${theta} -I 2 1 1 -ma x 0 5 x -ej 2.0 2 1 -T -L > msout
scrm 2 ${rep} -t ${theta} -I 2 1 1 -ma x 0 5 x -ej 2.0 2 1 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

