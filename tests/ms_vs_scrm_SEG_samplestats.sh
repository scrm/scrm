#!/bin/bash

mkdir test-SEGsamplestats
cd test-SEGsamplestats
rm *pdf


rep=100000

## compare SEG data
COMPAREFILE=compareSEG
rm ${COMPAREFILE}

theta=10

source ../chisq_r.src

source ../ks_r.src

source ../process_sample_stats.src

#case 1 
echo "4_samples" > current_case
ms 4 ${rep} -t ${theta} | sample_stats > ms_stats
scrm 4 ${rep} -t ${theta} | sample_stats > scrm_stats
foo


#case 2
echo "5_samples" > current_case
ms 5 ${rep} -t ${theta} | sample_stats > ms_stats
scrm 5 ${rep} -t ${theta} | sample_stats > scrm_stats
foo

#case 3
#2 sub population, 6 samples from subpopulation 1, and 5 samples from subpopulation 2, with rate 10 from 1 to 2, and rate 5 from 2 to 1
prefix="2groups6sam5sam_mig_x_10_5_x"
echo ${prefix} > current_case
ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 10.0 5.0 x | sample_stats > ms_stats
scrm 4 ${rep} -t ${theta} -I 2 2 2 -ma x 10.0 5.0 x | sample_stats > scrm_stats
foo

#case 4
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate matrix on manual page 5
prefix="3groups10sam4sam1sam_mig_x123x456x"
echo ${prefix} > current_case
ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x | sample_stats > ms_stats
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x | sample_stats > scrm_stats
foo
