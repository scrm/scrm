#!/bin/bash

dir=test-SUBPOP
mkdir ${dir}
cd ${dir}
rm *pdf


rep=100000

## compare population sturture for a single population data
COMPAREFILE=compareSubPop
rm ${COMPAREFILE}

theta=10

source ../chisq_r.src

source ../ks_r.src

source ../tmrca_r.src

source ../bl_r.src

source ../process_sample_stats.src	


#case 1 
echo "10_samples_case1" > current_case
rm ms* scrm*
ms 10 ${rep} -t ${theta} -I 2 2 8 -eN 0.4 10.01 -eN 1 0.01 -en 0.25 2 0.2 -ej 3 2 1 -T -L > msout
scrm 10 ${rep} -t ${theta} -I 2 2 8 -eN 0.4 10.01 -eN 1 0.01 -en 0.25 2 0.2 -ej 3 2 1 -T -L > scrmout


cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

#case 4
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate 5

echo "3groups10sam4sam1sam_mig_offdiag5" > current_case
rm ms* scrm*

	
ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -eN 0.8 15 -ej .7 2 1 -ej 1 3 1 -T -L > msout
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -eN 0.8 15 -ej .7 2 1 -ej 1 3 1 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo



#case 5
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate matrix on manual page 5

echo "3groups10sam4sam1sam_mig_x123x456x" > current_case

rm ms* scrm*
	
ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -eN 1 .1 -eN 3 10 -ej .7 2 1 -ej 4 3 1 -T -L > msout
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -eN 1 .1 -eN 3 10 -ej .7 2 1 -ej 4 3 1 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo
