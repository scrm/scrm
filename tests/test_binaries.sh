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

echo "Testing Initial Tree..."
debug_call 5 1 || exit 1
debug_call 3 100 -t 5 || exit 2 
debug_call 100 1 || exit 3 

echo "Testing Recombinations..."
debug_call 4 10 -r 5 100 -T || exit 4 
debug_call 6 10 -r 1 100 -t 5 || exit 5 

echo "Testing Pruning..."
debug_call 10 10 -r 100 1000 -l 10 || exit 7 
debug_call 3 10 -r 100 1000 -l 0 || exit 6 
