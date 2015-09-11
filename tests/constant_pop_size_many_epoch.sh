#!/bin/bash

function test_differ(){
	echo -n "."
	
    ./scrm $@ -seed ${seed} | tail -n +2 > tmp
	./scrm $@ -seed ${seed} -eN 0.000000 1 -eN 0.006475 1 -eN 0.041921 1 -eN 0.053825 1 -eN 0.067271 1 -eN 0.082459 1 -eN 0.099613 1 -eN 0.118989 1 -eN 0.140874 1 -eN 0.165593 1 -eN 0.193514 1 -eN 0.225051 1 -eN 0.260672 1 -eN 0.300906 1 -eN 0.346350 1 -eN 0.397680 1 -eN 0.455658 1 -eN 0.521144 1 -eN 0.595111 1 -eN 0.678657 1 -eN 0.773023 1 -eN 0.879609 1 -eN 1.000000 1 | tail -n +2 > tmpepoch1
	diff tmp tmpepoch1
	if [ $? -ne 0 ]; then
      echo ""
      echo "Failed at sample ${sample_size}, random seed:" ${seed}
	  echo "Check of First"
      exit 1
    fi	
	
	./scrm $@ -seed ${seed} | tail -n +2 > tmp
	./scrm $@ -seed ${seed} -eN 0.000000 1 -eN 0.006475 1 -eN 0.013789 1 -eN 0.022050 1 -eN 0.031381 1 -eN 0.041921 1 -eN 0.053825 1 -eN 0.067271 1 -eN 0.082459 1 -eN 0.099613 1 -eN 0.118989 1 -eN 0.140874 1 -eN 0.165593 1 -eN 0.193514 1 -eN 0.225051 1 -eN 0.260672 1 -eN 0.300906 1 -eN 0.346350 1 -eN 0.397680 1 -eN 0.455658 1 -eN 0.521144 1 -eN 0.595111 1 -eN 0.678657 1 -eN 0.773023 1 -eN 0.879609 1 -eN 1.000000 1 | tail -n +2 > tmpepoch1
	diff tmp tmpepoch1
	if [ $? -ne 0 ]; then
      echo ""
      echo "Failed at sample ${sample_size}, random seed:" ${seed}
	  echo "Check of Second"
      exit 1
    fi	
}

for sample_size in $( seq 2 30 ); do
	for seed in $( seq 1 100 ); do 
		test_differ ${sample_size} 1 -r 130 1000000 -l 0 -T || exit 1
		#../scrm sample_size ${sample_size} 1 -r 130 1000000 -seed ${seed} > 
	done
done
