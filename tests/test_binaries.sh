#!/bin/bash
#
# Test the binaries using scrms build in debug checks 
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Licence:  GPLv3 or later
#

function test_scrm {
  echo -n " scrm $@ "
  for i in `seq 1 10`; do
    echo -n "."

    # Test using scrm self-checks
    ./scrm_dbg $@ -seed $i > /dev/null 
    if [ $? -ne 0 ]; then
      echo ""
      echo "Executing \"./scrm_dbg $@ -seed $i\" failed."
      echo "Debug Call: make -mj2 scrm_dbg && ./scrm_dbg $@ -seed $i 2>&1 | less"
      exit 1
    fi

    # Test for memory leaks
    valgrind --error-exitcode=1 --leak-check=full -q ./scrm $@ -seed $i > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo "Valgrind check of \"./scrm $@ -seed $i\" failed."
      exit 1
    fi

  done
  echo " done."
}


echo "Testing Initial Tree"
 test_scrm 5 1000 -t 5 || exit 1
 test_scrm 3 1000 -t 5 -L -oSFS || exit 1
 test_scrm 100 20 -O || exit 1
 test_scrm 250 1 -T || exit 1
echo ""

echo "Testing Recombinations"
 test_scrm 4 500 -r 5 100 || exit 1
 test_scrm 6 500 -r 1 100 -t 5 -L -T || exit 1
 test_scrm 8 500 -r 1 100 -t 5 -oSFS -O || exit 1
echo ""

echo "Testing Pruning"
 test_scrm 10 200 -r 10 500 -l 10 -t 5 -oSFS || exit 1
 test_scrm 3 500 -r 10 500 -l 0 -t 1 -oSFS || exit 1
 test_scrm 10 200 -r 10 500 -l 2r -L || exit 1
 test_scrm 10 200 -r 1 500 -l -1 -T || exit 1
echo ""

echo "Testing Migration"
 test_scrm 5 50 -r 5 100 -I 2 3 2 1.2 || exit 1
 test_scrm 10 1 -r 20 200 -I 5 2 2 2 2 2 0.75 -l 5 || exit 1
 test_scrm 10 1 -r 10 100 -I 2 7 3 0.5 -eM 0.3 1.1 --print-model || exit 1
 test_scrm 10 1 -r 10 100 -I 2 7 3 -m 1 2 0.3 -em 0.5 2 1 0.6 -eM 2.0 1 || exit 1
 test_scrm 20 100 -I 3 2 2 2 1.0 -eI 1.0 2 2 2 -eI 2.0 2 3 3 || exit 1
echo ""

echo "Testing Size Change"
  test_scrm 10 2 -r 1 100 -I 3 3 3 4 0.5 -eN 0.1 0.05 -eN 0.2 0.5 --print-model || exit 1 
  test_scrm 10 2 -r 10 100 -I 3 3 3 4 0.5 -eN 0.1 0.05 -eN 0.2 0.5 -l 10 || exit 1 
echo ""

echo "Testing Splits & Merges"
 test_scrm 5 30 -r 20 200 -I 2 3 2 0.4 -ej 1.1 2 1 -l 25  || exit 1
 test_scrm 6 30 -r 20 200 -I 3 2 2 2 -ej 0.2 2 1 -ej 0.25 3 1 -l 25 || exit 1
 test_scrm 20 2 -r 5 200 -I 2 10 10 1.5 -es 1.6 2 0.5 -eM 2.0 1 -l 25 || exit 1
 test_scrm 20 2 -r 5 200 -I 2 10 10 1.5 -es 0.9 1 0.8 -es 1.6 2 0.5 -eM 2.0 1 -l 50 || exit 1
 test_scrm 20 2 -r 5 200 -I 2 10 10 -eps 1.0 2 1 0.0 || exit 1
echo ""

echo "Testing Growth"
 test_scrm 5 100 -r 20 200 -G 1.5 -l 25 || exit 1
 test_scrm 5 60 -r 2 200 -G -2.5 -eG 1 0.0 -l 25 || exit 1
 test_scrm 4 30 -r 2 200 -G -2.5 -eN 1 0.25 -eG 2 0.0 -eN 2.5 0.25 -l 25 || exit 1
echo ""

echo "Testing Variable Rates"
 test_scrm 3 200 -r 2 100 -t 5 -st 10 10 -sr 20 5 -st 30 1 -sr 40 0 -st 50 20 -T || exit 1
echo ""

echo "Various Edge Cases"
 test_scrm 6 1 -I 2 3 3 0.5 -r 1 100 -es 0 2 0.5 -ej 1 3 1 -t 1  || exit 1
 test_scrm 10 10 -es 1.0 1 0.5 -ej 1.0 2 1 || exit 1
echo ""

