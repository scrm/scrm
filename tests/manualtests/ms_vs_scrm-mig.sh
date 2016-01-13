#!/bin/bash

#Compare summary statistics of ms and scrm for migration

mkdir test-mig
cd test-mig
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
#2 sub population, 2 samples from each subpopulation, mutation rate is 5
echo "2groups2sam2sam_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 4 ${rep} -t ${theta} -I 2 2 2 5.0 -T -L > msout
../../../scrm 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

#case 2
#2 sub population, 6 samples from subpopulation 1, and 5 samples from subpopulation 2, with rate 10 from 1 to 2, and rate 5 from 2 to 1

echo "2groups6sam5sam_mig_x_10_5_x" > current_case
rm ms* scrm*

ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 10.0 5.0 x -T -L > msout
scrm 4 ${rep} -t ${theta} -I 2 2 2 -ma x 10.0 5.0 x -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo


#case 3
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate 5

echo "3groups10sam4sam1sam_mig_sym10" > current_case
rm ms* scrm*

ms 15 ${rep} -t ${theta} -I 3 10 4 1 10.0 -T -L > msout
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

#case 4
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate 5

echo "3groups10sam4sam1sam_mig_offdiag5" > current_case
rm ms* scrm*


ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -T -L > msout
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo



#case 5
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate matrix on manual page 5

echo "3groups10sam4sam1sam_mig_x123x456x" > current_case

rm ms* scrm*

ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -T -L > msout
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo
