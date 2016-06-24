#/bin/bash
wget https://webshare.uchicago.edu/users/rhudson1/Public/ms.folder/ms.tar.gz -o /dev/null
tar -xf ms.tar.gz
cd msdir; gcc -o ms ms.c streec.c rand1.c -lm ; gcc -o sample_stats sample_stats.c tajd.c -lm; cd ..

wget --no-check-certificate https://github.com/hybridLambda/hybrid-Lambda/archive/release.tar.gz -o /dev/null
tar -xf release.tar.gz
cd hybrid-Lambda-release ; ./bootstrap ; make hybrid-Lambda ; cd ..

make scrm

mkdir test-demo
cd test-demo
rm *pdf

rep=1000000
recombRep=10000
seqlen=100000
theta=10
r=10


testingFunction() {
  echo "testing $1 $2 $3 $4"
  R --slave "--args $1 $2 $3 $4" < ../tests/manualtests/ksTest.r > test.out
  grep "Status Ok" test.out
  if [ $? -ne 0 ]; then
    echo ""
    echo "Failed testing for $1"
    exit 1
  fi
}

computeMoments(){
  program=$1
  cat ${program}out | gawk '/^\/\//{f="xx"++d} f{print > f} '
  for file in $(seq 1 1 ${recombRep})
    do
    grep ";" xx${file} | sed -e 's/\[.*\]//g' > xxTrees
    ../hybrid-Lambda-release/hybrid-Lambda -gt xxTrees -bl -tmrca -o xx${file}Trees 2> /dev/null
    grep ";" xx${file} | sed -e 's/\[//g' -e 's/\].*;//g' > xx${file}Trees_freq
    done
  R --slave "--args ${recombRep} ${seqlen} $1" < ../tests/manualtests/tmrcaMoments.r > tmrcaMoments.out
  find . -name "xx*" -print0 | xargs -0 rm
}

tmrca_no_recomb(){
  cat msout | grep "time:" > mstime
  cat scrmout | grep "time:" >  scrmtime

  cut -f 2 mstime > mstmrca
  cut -f 2 scrmtime > scrmtmrca
  testingFunction "TMRCA" $1 "mstmrca" "scrmtmrca" || exit 1

  cut -f 3 mstime > msbl
  cut -f 3 scrmtime > scrmbl
  testingFunction "BL" $1 "msbl" "scrmbl" || exit 1
}

moment_recomb(){
  computeMoments ms
  computeMoments scrm

  objs=("tmrca" "bl")
  for obj in ${objs[@]}
    do
    for moment in  $(seq 1 1 4)
      do
      cut -f ${moment} ms${obj}_moments > msmoment
      cut -f ${moment} scrm${obj}_moments > scrmmoment
      testingFunction ${obj}${moment} $1 "msmoment" "scrmmoment" || exit 1
      done
    done
}


testingScript(){
  cat msout | ../msdir/sample_stats > ms_stats
  cat scrmout | ../msdir/sample_stats > scrm_stats

  cut -f 6 ms_stats > msdata
  cut -f 6 scrm_stats > scrmdata
  testingFunction "Tajima_D" $1 "msdata" "scrmdata" || exit 1

  cut -f 2 ms_stats > msdata
  cut -f 2 scrm_stats > scrmdata
  testingFunction "Pairwise_difference" $1 "msdata" "scrmdata" || exit 1

  cut -f 8 ms_stats > msdata
  cut -f 8 scrm_stats > scrmdata
  testingFunction "theta_H" $1 "msdata" "scrmdata" || exit 1

  cut -f 10 ms_stats > msdata
  cut -f 10 scrm_stats > scrmdata
  testingFunction "H" $1 "msdata" "scrmdata" || exit 1
}


runTest(){
	label="$1"
	secondTestingScript="$2"
	commandLine="$3"
	# run ms and scrm to produce the output
	../msdir/ms ${commandLine} -seed 1 1 1 > msout
	../scrm ${commandLine} -seed 1 1 1 > scrmout
	# run test
	succeed=1
	testingScript "${label}" || succeed=0
	${secondTestingScript} "${label}" || succeed=0
	if [ "${succeed}"="0" ] ; then
		echo "Failure on first round -- trying again with different seed"
		../msdir/ms ${commandLine} -seed 2 2 2 > msout
		../scrm ${commandLine} -seed 2 2 2 > scrmout
		testingScript "${label}" || exit 1
		${secondTestingScript} "${label}" || exit 1
		echo "Success on second round"
	fi
}


#case 1
#2 sub population, 2 samples from each subpopulation, mutation rate is 5
rm ms* scrm*

runTest "2groups2sam2sam_mig5" tmrca_no_recomb "4 ${rep} -t ${theta} -I 2 2 2 5.0 -T -L"
runTest "2groups3sam3sam_mig5" tmrca_no_recomb "6 ${rep} -t ${theta} -I 2 3 3 5.0 -T -L"
runTest "8sam_recomb" moment_recomb "8 ${recombRep} -t ${theta} -r ${r} ${seqlen} -T"
runTest "10sam_recomb" moment_recomb "10 ${recombRep} -t ${theta} -r ${r} ${seqlen} -T"


# Cleaning up
cd ..

rm -r ms.tar.gz msdir
rm -r release.tar.gz hybrid-Lambda-release
rm -r test-demo
