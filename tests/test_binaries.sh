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
    ./src/scrm_dbg $@ -seed $i > /dev/null 
    if [ $? -ne 0 ]; then
      echo "Executing \"./src/scrm_dbg $@ -seed $i\" failed."
      exit 1
    fi
  done
}

function normal_call {
  for i in `seq 1 100`; do
    loci=`./src/scrm $@ -seed $i | grep -c "//"`
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
debug_call 3 100 -t 5 || exit 2 
debug_call 100 1 -T || exit 3 

echo "Testing Recombinations..."
debug_call 4 10 -r 5 100 -T || exit 4 
debug_call 6 10 -r 1 100 -t 5 || exit 5 

echo "Testing Pruning..."
debug_call 10 10 -r 10 500 -l 10 || exit 6
debug_call 3 10 -r 10 500 -l 0 || exit 7

echo "Testing Migration..."
debug_call 10 10 -r 50 500 -I 2 5 5 0.5 -l 100 || exit 8
debug_call 10 10 -r 50 500 -I 2 5 5 -M 1.7 -l 100  || exit 9

echo "Testing Split..."
debug_call 10 10 -r 50 500 -I 2 5 5 0.5 -ej 1.0 2 1 -l 100 || exit 10
