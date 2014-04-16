#!/bin/bash
#
# Test the binaries using scrms build in debug checks 
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-10-10
# Licence:  GPLv3 or later
#

function debug_call {
  for i in `seq 1 100`; do
    ./scrm_dbg $@ -seed $i > /dev/null 
    if [ $? -ne 0 ]; then
      echo "Executing \"./src/scrm_dbg $@ -seed $i\" failed."
      exit 1
    fi
  done
}

function normal_call {
  for i in `seq 1 100`; do
    loci=`./scrm $@ -seed $i | grep -c "//"`
    if [ $loci -ne $2 ]; then
      echo "Executing \"./src/scrm $@ -seed $i\" failed."
      exit 1
    fi
  done
}

echo "Testing Binary..."
normal_call 8 1 -t 5 || exit 1
normal_call 9 10 -t 5 -T || exit 1
normal_call 10 11 -r 10 100 -t 5 || exit 1

echo "Testing Initial Tree..."
debug_call 5 1 -t 5 || exit 1
debug_call 3 100 -t 5 || exit 1
debug_call 100 1 -T || exit 1

echo "Testing Recombinations..."
debug_call 4 10 -r 5 100 -T || exit 1
debug_call 6 10 -r 1 100 -t 5 || exit 1

echo "Testing Pruning..."
debug_call 10 10 -r 10 500 -l 10 -t 5 -oSFS || exit 1
debug_call 3 10 -r 10 500 -l 0 -t 5 -oSFS || exit 1

echo "Testing Migration..."
debug_call 5 2 -r 5 100 -I 2 3 2 1.2 || exit 1
debug_call 5 2 -r 20 200 -I 2 3 2 1.2 -l 25 || exit 1

echo "Testing Split..."
debug_call 5 5 -r 20 200 -I 2 3 2 0.4 -ej 1.1 2 1 -l 25  || exit 1
debug_call 6 2 -r 20 200 -I 3 2 2 2 -ej 0.2 2 1 -ej 0.25 3 1 -l 25 || exit 1

echo "Testing Growth..."
debug_call 5 20 -r 20 200 -G 0.5 -l 25 || exit 1
debug_call 5 10 -r 2 200 -G -2.5 -eG 1 0.0 -l 25  || exit 1

echo "Testing Variable Rates..."
debug_call 2 1 -r 2 100 -t 5 -st 10 10 -sr 20 5 -st 30 1 -sr 40 0 -st 50 20 -T || exit 1
