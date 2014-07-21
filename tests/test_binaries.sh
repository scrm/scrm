#!/bin/bash
#
# Test the binaries using scrms build in debug checks 
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-10-10
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
      echo "Debug Call: ./scrm_dbg $@ -seed $i 2>$1 | less"
      exit 1
    fi

    # Test for memory leaks
    valgrind --error-exitcode=1 --leak-check=full -q ./scrm $@ -seed $i > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo "Valgrind check of \"./scrm $@ -seed $i\" failed."
      exit 1
    fi

    # Test for reproducibility
    ./scrm $@ -seed $i
    ./scrm $@ -seed $i
    hash_1=$(./scrm $@ -seed $i | sha512sum)
    hash_2=$(./scrm $@ -seed $i | sha512sum)
    if [ "$hash_1" != "$hash_2" ]; then
      echo ""
      echo "Not reproducible: \"./scrm $@ -seed $i\"" 
      echo "hash1: $hash_1"
      echo "hash2: $hash_2"
      exit 1
    fi
  done
  echo " done."
}

echo "Testing Initial Tree"
 test_scrm 5 1000 -t 5 || exit 1
 test_scrm 3 1000 -t 5 -L -oSFS || exit 1
 test_scrm 100 20 -T || exit 1
 test_scrm 250 1 || exit 1
echo ""

echo "Testing Recombinations"
 test_scrm 4 500 -r 5 100 || exit 1
 test_scrm 6 500 -r 1 100 -t 5 -L || exit 1
 test_scrm 8 500 -r 1 100 -t 5 -oSFS || exit 1
echo ""

echo "Testing Pruning"
 test_scrm 10 200 -r 10 500 -l 10 -t 5 -oSFS || exit 1
 test_scrm 3 500 -r 10 500 -l 0 -t 5 -oSFS || exit 1
echo ""

echo "Testing Migration"
 test_scrm 5 50 -r 5 100 -I 2 3 2 1.2 || exit 1
 test_scrm 10 1 -r 20 200 -I 5 2 2 2 2 2 0.75 -l 5 || exit 1
 test_scrm 20 10 -r 5 100 -I 3 2 2 2 1.0 -eI 1.0 2 2 2 -eI 2.0 2 3 3 || exit 1
echo ""

echo "Testing Split"
 test_scrm 5 30 -r 20 200 -I 2 3 2 0.4 -ej 1.1 2 1 -l 25  || exit 1
 test_scrm 6 30 -r 20 200 -I 3 2 2 2 -ej 0.2 2 1 -ej 0.25 3 1 -l 25 || exit 1
echo ""

echo "Testing Growth..."
 test_scrm 5 100 -r 20 200 -G 1.5 -l 25 || exit 1
 test_scrm 5 60 -r 2 200 -G -2.5 -eG 1 0.0 -l 25  || exit 1
echo ""

echo "Testing Variable Rates..."
 test_scrm 3 200 -r 2 100 -t 5 -st 10 10 -sr 20 5 -st 30 1 -sr 40 0 -st 50 20 -T || exit 1
echo ""
